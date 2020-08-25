;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************


undef("regrid_data")
function regrid_data(var,interp_method,lat2d,lon2d,lat,lon,src_mask,dst_mask)
local var_regrid
begin

    Opt                   = True

    Opt@SrcFileName       = "WRF_SCRIP.nc"      ; Name of source and ???
    Opt@DstFileName       = "1degree.nc"    ; destination files  ???

    Opt@InterpMethod      = interp_method

    Opt@SrcRegional       = True
    Opt@SrcGridMask       = src_mask

    Opt@DstRegional       = True
    ; Opt@DstLLCorner       = (/ minlat,minlon /)
    ; Opt@DstURCorner       = (/ maxlat,maxlon /)
    Opt@DstGridType       = "rectilinear"
    Opt@DstGridLat        = lat
    Opt@DstGridLon        = lon
    Opt@DstGridMask       = dst_mask

    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@Debug             = True

    var@lat2d = lat2d              ; This information will be used by
    var@lon2d = lon2d              ; ESMF_regrid for the source grid
    var_regrid = ESMF_regrid(var,Opt)     ; Do the regridding for TMP
    printVarSummary(var_regrid)
    return(var_regrid)
end

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; =================================== Path =====================================
path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/Princeton_ctl/LIS_output"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl_para/outputs"
path3 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl/outputs"


; ================================= Predefine ==================================
year_s = 2003
year_e = 2008
year_sum = year_e-year_s+1

Soil_thickness1 =(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE
betaroot       = (/ \
                  0.943, 0.962, 0.966, 0.961, 0.966, \
                  0.914, 0.964, 0.972, 0.943, 0.943, \
                  0.961, 0.961, 0.961, 0.961, 0.961, \
                  0.975, 0.975 /)

; ================================ Set Regrid ==================================
interp_method = "neareststod"  ;"neareststod"; "conserve" "patch";"bilinear"

filename  = path1+"/LIS.CABLE."+year_s+"090100.d01.nc"
f = addfile (filename,"r")
filename1 = path2+"/cable_out_"+year_s+".nc"
f1 = addfile (filename1,"r")


lat2d    = f->lat(0,:,:)
lon2d    = f->lon(0,:,:)
lon2d    = where(lon2d.lt.0,360+lon2d,lon2d)

minlon  = ceil(min(lon2d))+0.5 ; 116.5516
minlat  = ceil(min(lat2d))+0.5 ; -44.43771
maxlon  = floor(max(lon2d))-0.5 ; 154.5981
maxlat  = floor(max(lat2d))-0.5 ; -12.5853

src_mask = where(.not.ismissing(f->Landcover_inst(0,:,:)),1,0)
dst_mask = where(.not.ismissing(f1->Qs(0,{minlat:maxlat},{minlon:maxlon})),1,0)

lat      = fspan(minlat,maxlat,dimsizes(dst_mask(:,0)))
lon      = fspan(minlon,maxlon,dimsizes(dst_mask(0,:)))
delete(f)
delete(f1)

; study region
A = minlat
B = maxlat
C = minlon
D = maxlon


; ================================ READ CABLE ==================================
filename  = path1+"/LIS.CABLE."+year_s+"010100.d01.nc"
f         = addfile (filename,"r")
Landcover = f->Landcover_inst(0,:,:)
RelSMC    = dim_avg_n_Wrap(f->RelSMC_tavg,0)
Fwsoil_LIS= dim_avg_n_Wrap(f->fwsoil_tavg,0)
RelSMC    = 0.
Fwsoil_LIS= 0.
delete(f)

filename1 = path2+"/cable_out_"+year_s+".nc"
f1        = addfile (filename1,"r")
Fwsoil_HESS = dim_sum_n_Wrap(f1->Fwsoil(:,{A:B},{C:D}),0)
Fwsoil_HESS = 0.
delete(f1)

filename2 = path3+"/cable_out_"+year_s+".nc"
f2        = addfile (filename2,"r")
Fwsoil_HESS_def = dim_sum_n_Wrap(f2->Fwsoil(:,{A:B},{C:D}),0)
Fwsoil_HESS_def = 0.
delete(f2)

do year = year_s,year_e
    print("year = " +year)
    ; _________________________ read LIS-CABLE ____________________________
    do month = 1, 12
      if month .lt. 10 then
        filename  = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename  = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      f = addfile (filename,"r")
      RelSMC     = RelSMC   + dim_avg_n_Wrap(f->RelSMC_tavg,0)
      Fwsoil_LIS = Fwsoil_LIS + dim_avg_n_Wrap(f->fwsoil_tavg,0)
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
Fwsoil_LIS   = Fwsoil_LIS/year_sum/12
Fwsoil_HESS  = Fwsoil_HESS/year_sum
Fwsoil_HESS_def  = Fwsoil_HESS_def/year_sum

; ; _________________________ Calculate Froot ________________________________
; froot = new((/17,6/),float)
; froot = 0.
;
; do vegtype = 0, 16
;   depth = 0.
;   do soil = 0,5
;      print(betaroot(vegtype))
;      depth = depth + (Soil_thickness1(soil)*100.0)
;      print(depth)
;      print(froot(vegtype,soil))
;      froot(vegtype,soil) = 1. - betaroot(vegtype)^depth
;   end do
;   do soil = 5,1,1
;      froot(vegtype,soil) = froot(vegtype,soil) - froot(vegtype,soil-1)
;   end do
;   print(froot(vegtype,:))
; end do
;
; Froot_LIS = RelSMC
; Froot_LIS = 0.
; printVarSummary(Froot_LIS)
;
; do soil = 0, 5
;   do vegtype = 1, 17
;     Froot_LIS(soil, :, :) = where(Landcover .eq. vegtype, froot(vegtype-1 ,soil), Froot_LIS(soil, :, :))
;   end do
; end do
;
; Fwsoil_LIS_cal = dim_sum_n_Wrap(Froot_LIS * RelSMC, 0)
; printVarSummary(Fwsoil_LIS_cal)

; ================================== Regrid ====================================

; Fwsoil_LIS_cal_regrid = regrid_data(Fwsoil_LIS_cal,interp_method, lat2d, lon2d, \
;                         lat, lon, src_mask, dst_mask )
Fwsoil_LIS_regrid     = regrid_data(Fwsoil_LIS,interp_method, lat2d, lon2d, \
                        lat, lon, src_mask, dst_mask )


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
gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")             ;����ɫ��

; Fwsoil
plots = new(5,graphic)
pres@txString = "Fwsoil (-)"
levels = fspan(0,1,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,Fwsoil_LIS_regrid,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,Fwsoil_HESS,res)
res@tiMainString = "CABLE(HESS) def param"
plots(2) = gsn_csm_contour(wks1,Fwsoil_HESS_def,res)

levels = fspan(-0.5,0.5,21)
res@cnLevels = levels
res@cnFillPalette = "ViBlGrWhYeOrRe"
res@tiMainString = "LIS-CABLE - CABLE(HESS)"
plots(3) = gsn_csm_contour(wks1,Fwsoil_LIS_regrid-Fwsoil_HESS,res)
res@tiMainString = "CABLE(HESS) def param - CABLE(HESS)"
plots(4) = gsn_csm_contour(wks1,Fwsoil_HESS_def-Fwsoil_HESS,res)
; res@tiMainString = "LIS-CABLE_cal"
; plots(1) = gsn_csm_contour(wks1,Fwsoil_LIS_cal_regrid,res)
gsn_panel(wks1,(/plots/),(/2,3/),pres)

end
