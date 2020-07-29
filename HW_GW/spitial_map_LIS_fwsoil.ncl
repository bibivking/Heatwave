;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; =================================== Path =====================================
path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/nuwrf_25km/LISWRF_configs/LIS_output"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl/outputs"
path3 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/litter/outputs"

; ================================= Predefine ==================================
year_s = 2003
year_e = 2008
year_sum = year_e-year_s+1

; study region
A = -43.4377
B = -13.5853
C = 117.552
D = 153.598

Soil_thickness1 =(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE
betaroot       = (/ \
                  0.943, 0.962, 0.966, 0.961, 0.966, \
                  0.914, 0.964, 0.972, 0.943, 0.943, \
                  0.961, 0.961, 0.961, 0.961, 0.961, \
                  0.975, 0.975 /)


; ================================ READ CABLE ==================================
filename = path1+"/LIS.CABLE."+year_s+"010100.d01.nc"
f        = addfile (filename,"r")

Landcover = f->Landcover_inst(0,:,:)
RelSMC    = dim_avg_n_Wrap(f->RelSMC_inst,0)
RelSMC    = 0.
printVarSummary(RelSMC)
delete(f)

filename1 = path2+"/cable_out_"+year_s+".nc"
f1        = addfile (filename1,"r")

Fwsoil_HESS = dim_sum_n_Wrap(f1->Fwsoil(:,{A:B},{C:D}),0)
Fwsoil_HESS = 0.
printVarSummary(Fwsoil_HESS)
delete(f1)

filename2 = path3+"/cable_out_"+year_s+".nc"
f2        = addfile (filename2,"r")

Fwsoil_HESS_def = dim_sum_n_Wrap(f2->Fwsoil(:,{A:B},{C:D}),0)
Fwsoil_HESS_def = 0.
printVarSummary(Fwsoil_HESS_def)
delete(f2)

do year = year_s,year_e
    print("year = " +year)
    ; _________________________ read LIS-CABLE ____________________________
    do month = 1, 12
      if month .lt. 10 then
        filename = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      f = addfile (filename,"r")
      RelSMC     = RelSMC   + dim_avg_n_Wrap(f->RelSMC_inst,0)
      delete(f)
    end do ; month

    ; _________________________ read CABLE (HESS) _________________________
    filename1 = path2+"/cable_out_"+year+".nc"
    f1 = addfile (filename1,"r")
    Fwsoil_HESS = Fwsoil_HESS   + dim_avg_n_Wrap(f1->Fwsoil(:,{A:B},{C:D}),0)
    delete(f1)

    ; _________________________ read CABLE (HESS) _________________________
    filename2 = path3+"/cable_out_"+year+".nc"
    f2 = addfile (filename2,"r")
    Fwsoil_HESS_def = Fwsoil_HESS_def   + dim_avg_n_Wrap(f2->Fwsoil(:,{A:B},{C:D}),0)
    delete(f2)
end do ; year

RelSMC       = RelSMC/year_sum/12
Fwsoil_HESS  = Fwsoil_HESS/year_sum
Fwsoil_HESS_def  = Fwsoil_HESS_def/year_sum

; _________________________ Calculate Froot ________________________________
froot = new((/17,6/),float)
froot = 0.

do vegtype = 0, 16
  depth = 0.
  do soil = 0,5
     print(betaroot(vegtype))
     depth = depth + (Soil_thickness1(soil)*100.0)
     print(depth)
     print(froot(vegtype,soil))
     froot(vegtype,soil) = 1. - betaroot(vegtype)^depth
  end do
  do soil = 5,1,1
     froot(vegtype,soil) = froot(vegtype,soil) - froot(vegtype,soil-1)
  end do
  print(froot(vegtype,:))
end do

Froot_LIS = RelSMC
Froot_LIS = 0.
printVarSummary(Froot_LIS)

do soil = 0, 5
  do vegtype = 1, 17
    Froot_LIS(soil, :, :) = where(Landcover .eq. vegtype, froot(vegtype-1 ,soil), Froot_LIS(soil, :, :))
  end do
end do

Fwsoil_LIS = dim_sum_n_Wrap(Froot_LIS * RelSMC, 0)
printVarSummary(Fwsoil_LIS)


; ================================= Plotting ===================================
res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0

res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ
res@lbLabelBarOn       = True

;************** ����labelbar ***************
res@lbBoxLinesOn       = True                       ;�ر�lbar box ����
res@lbTitleFont        = 25
res@lbLabelFont        = 25
res@lbTitleFontHeightF = 0.013
res@lbLabelFontHeightF = 0.013
res@txString      = ""
res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont      = 25

;*************** ���õ�ֵ�� *****************
res@cnLinesOn          = False                       ; �رյ�ֵ������
res@cnLineColor        = "black"
res@cnLineThicknessF   = 1.5
res@cnLineLabelsOn     = False
res@gsnMaximize        = True
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode = "ExplicitLevels"

;************ ����������ֵ��Χ **************
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015
res@tmYLMode  = "Explicit"
res@tmXBMode  = "Explicit"

;**************** ͼ�������� *****************
pres                    = True                                      ; ���� panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
;pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01

res@gsnLeftString = ""
res@gsnRightString = ""


; ============================= plot water balance =============================
pic1 = "./plots/LIS_fwsoil"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

; Fwsoil
plots = new(3,graphic)
pres@txString = "Fwsoil (-)"
levels = fspan(0,1,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,Fwsoil_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,Fwsoil_HESS,res)
res@tiMainString = "CABLE(HESS)-litter"
plots(2) = gsn_csm_contour(wks1,Fwsoil_HESS_def,res)
gsn_panel(wks1,(/plots/),(/1,3/),pres)

end
