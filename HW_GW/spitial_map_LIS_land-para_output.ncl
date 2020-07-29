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
path3 = "/g/data/w35/mm3972/model/cable/runs/my_version/outputs/GSWP"

; ================================= Predefine ==================================
year_s = 2003
year_e = 2008
year_sum = year_e-year_s+1

; study region
A = -43.4377
B = -13.5853
C = 117.552
D = 153.598

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Soil_thickness1 =(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE

; ================================ READ CABLE ==================================
filename = path1+"/LIS.CABLE."+year_s+"010100.d01.nc"
f        = addfile (filename,"r")

Albedo_LIS     = dim_avg_n_Wrap(f->Albedo_inst,0)
Landcover_LIS  = dim_avg_n_Wrap(f->Landcover_inst,0)
Soiltype_LIS   = dim_avg_n_Wrap(f->Soiltype_inst,0)
sand_LIS       = dim_avg_n_Wrap(f->SandFrac_inst,0)
clay_LIS       = dim_avg_n_Wrap(f->ClayFrac_inst,0)
silt_LIS       = dim_avg_n_Wrap(f->SiltFrac_inst,0)

filename1 = path2+"/cable_out_"+year_s+".nc"
f1        = addfile (filename1,"r")

Albedo_HESS = dim_avg_n_Wrap(f1->Albedo(:,{A:B},{C:D}),0)
iveg_HESS   = f1->iveg({A:B},{C:D})
isoil_HESS  = f1->isoil({A:B},{C:D})
sand_HESS   = (f1->sand(0,{A:B},{C:D})*Soil_thickness(0) +\
              f1->sand(1,{A:B},{C:D})*Soil_thickness(1) +\
              f1->sand(2,{A:B},{C:D})*Soil_thickness(2) +\
              f1->sand(3,{A:B},{C:D})*Soil_thickness(3) +\
              f1->sand(4,{A:B},{C:D})*Soil_thickness(4) +\
              f1->sand(5,{A:B},{C:D})*Soil_thickness(5))/sum(Soil_thickness)
silt_HESS   = (f1->silt(0,{A:B},{C:D})*Soil_thickness(0) +\
              f1->silt(1,{A:B},{C:D})*Soil_thickness(1) +\
              f1->silt(2,{A:B},{C:D})*Soil_thickness(2) +\
              f1->silt(3,{A:B},{C:D})*Soil_thickness(3) +\
              f1->silt(4,{A:B},{C:D})*Soil_thickness(4) +\
              f1->silt(5,{A:B},{C:D})*Soil_thickness(5))/sum(Soil_thickness)
clay_HESS   = (f1->clay(0,{A:B},{C:D})*Soil_thickness(0) +\
              f1->clay(1,{A:B},{C:D})*Soil_thickness(1) +\
              f1->clay(2,{A:B},{C:D})*Soil_thickness(2) +\
              f1->clay(3,{A:B},{C:D})*Soil_thickness(3) +\
              f1->clay(4,{A:B},{C:D})*Soil_thickness(4) +\
              f1->clay(5,{A:B},{C:D})*Soil_thickness(5))/sum(Soil_thickness)

path3 = "/g/data/w35/mm3972/model/cable/runs/my_version/outputs/GSWP"
filename2 = path3+"/cable_out_"+year_s+".nc"
f2        = addfile(filename2,"r")

Albedo_HESS_def = dim_avg_n_Wrap(f2->Albedo(:,{A:B},{C:D}),0)

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
pic1 = "./plots/LIS_land_para_output"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

; Albedo
plots = new(3,graphic)
pres@txString = "Albedo"
levels = fspan(0,0.4,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,Albedo_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,Albedo_HESS,res)
res@tiMainString = "CABLE(HESS)-def"
plots(2) = gsn_csm_contour(wks1,Albedo_HESS_def,res)
gsn_panel(wks1,(/plots/),(/1,3/),pres)
delete(plots)
delete(levels)
delete(res@cnLevels)

; Landcover
plots = new(2,graphic)
pres@txString = "Landcover (mm/y)"
levels = fspan(1,21,22)
res@cnLevels = levels
print(min(Landcover_LIS))
print(max(Landcover_LIS))
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,Landcover_LIS,res)
res@tiMainString = "CABLE(HESS)"
print(min(iveg_HESS))
print(max(iveg_HESS))
plots(1) = gsn_csm_contour(wks1,iveg_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

; Soil type
pres@txString = "Soil type"
levels = fspan(1,21,22)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
print(min(Soiltype_LIS))
print(max(Soiltype_LIS))
plots(0) = gsn_csm_contour(wks1,Soiltype_LIS,res)
res@tiMainString = "CABLE(HESS)"
print(min(isoil_HESS))
print(max(isoil_HESS))
plots(1) = gsn_csm_contour(wks1,isoil_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

delete(levels)
delete(res@cnLevels)

; Sand
pres@txString = "sand"
levels = fspan(0,1.,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,sand_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,sand_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

; Silt
pres@txString = "silt"
levels = fspan(0,1.,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,silt_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,silt_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

; Clay
pres@txString = "clay"
levels = fspan(0,1.,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,clay_LIS,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,clay_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

delete(plots)
delete(levels)
delete(res@cnLevels)

end
