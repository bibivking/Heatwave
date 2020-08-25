;*******************************************
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
path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/Princeton_ctl_para/LIS_output"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl_para/outputs"
path3 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl/outputs"
path4 = "/g/data/w35/mm3972/model/cable/src/CABLE-AUX/offline"


; ================================= Predefine ==================================
year_s = 2003
year_e = 2008
year_sum = year_e-year_s+1

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Soil_thickness1 = Soil_thickness
;(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE

; ================================ Set Regrid ==================================
interp_method = "neareststod"  ;"neareststod"; "conserve" "patch";"bilinear"

filename  = path1+"/LIS.CABLE."+year_e+"090100.d01.nc"
f = addfile (filename,"r")
filename1 = path2+"/cable_out_"+year_e+".nc"
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

; READ LIS-CABLE
filename = path1+"/LIS.CABLE."+year_e+"010100.d01.nc"
f        = addfile (filename,"r")

; Albedo_LIS    = dim_avg_n_Wrap(f->Albedo_tavg,0)
Landcover_LIS = dim_avg_n_Wrap(f->Landcover_inst,0)
Soiltype_LIS  = dim_avg_n_Wrap(f->Soiltype_inst,0)

sand_LIS      = dim_avg_n_Wrap(f->SandFrac_inst,0)
clay_LIS      = dim_avg_n_Wrap(f->ClayFrac_inst,0)
silt_LIS      = dim_avg_n_Wrap(f->SiltFrac_inst,0)
ssat_LIS      = dim_avg_n_Wrap(f->Porosity_inst,0)
bch_LIS       = dim_avg_n_Wrap(f->bch_inst,0)
sfc_LIS       = dim_avg_n_Wrap(f->sfc_inst,0)
swilt_LIS     = dim_avg_n_Wrap(f->swilt_inst,0)
hyds_LIS      = dim_avg_n_Wrap(f->hyds_inst,0)
sucs_LIS      = dim_avg_n_Wrap(f->sucs_inst,0)
css_LIS       = dim_avg_n_Wrap(f->css_inst,0)
rhosoil_LIS   = dim_avg_n_Wrap(f->rhosoil_inst,0)
froot_LIS     = ( dim_avg_n_Wrap(f->froot_inst(:,0,:,:),0)*Soil_thickness1(0) \
              + dim_avg_n_Wrap(f->froot_inst(:,1,:,:),0)*Soil_thickness1(1) \
              + dim_avg_n_Wrap(f->froot_inst(:,2,:,:),0)*Soil_thickness1(2) \
              + dim_avg_n_Wrap(f->froot_inst(:,3,:,:),0)*Soil_thickness1(3) \
              + dim_avg_n_Wrap(f->froot_inst(:,4,:,:),0)*Soil_thickness1(4) \
              + dim_avg_n_Wrap(f->froot_inst(:,5,:,:),0)*Soil_thickness1(5))/sum(Soil_thickness1)
elev_LIS      = dim_avg_n_Wrap(f->Elevation_inst,0)

; fwsoil_LIS    = dim_avg_n_Wrap(f->fwsoil_inst,0)
; LAI_LIS       = dim_avg_n_Wrap(f->LAI_inst,0)

; READ CABLE(HESS)
filename1 = path2+"/cable_out_"+year_s+".nc"
f1        = addfile (filename1,"r")

; Albedo_HESS = f1->Albedo(0,{A:B},{C:D}) ; Jan
iveg_HESS   = f1->iveg({A:B},{C:D})
isoil_HESS  = f1->isoil({A:B},{C:D})
elev_HESS   = f1->elev({A:B},{C:D})
; fwsoil_HESS = dim_avg_n_Wrap(f1->Fwsoil(:,{A:B},{C:D}),0)
; LAI_HESS    = dim_avg_n_Wrap(f1->LAI(:,{A:B},{C:D}),0)

sand_HESS   = f1->sand(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
silt_HESS   = f1->silt(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
clay_HESS   = f1->clay(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
ssat_HESS   = f1->ssat(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
bch_HESS    = f1->bch(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
sfc_HESS    = f1->sfc(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
swilt_HESS  = f1->swilt(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
hyds_HESS   = f1->hyds(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
sucs_HESS   = f1->sucs(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
css_HESS    = f1->css(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
rhosoil_HESS= f1->rhosoil(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)
froot_HESS  = f1->froot(0,{A:B},{C:D})*Soil_thickness(0)/sum(Soil_thickness)

; fwsoil_HESS = f1->Fwsoil(:,{A:B},{C:D})
; LAI_HESS    = f1->LAI(:,{A:B},{C:D})


do soil_layer = 1,5
  sand_HESS   = sand_HESS + f1->sand(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  silt_HESS   = silt_HESS + f1->silt(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  clay_HESS   = clay_HESS + f1->clay(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  ssat_HESS   = ssat_HESS + f1->ssat(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  bch_HESS    = bch_HESS + f1->bch(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  sfc_HESS    = sfc_HESS + f1->sfc(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  swilt_HESS  = swilt_HESS + f1->swilt(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  hyds_HESS   = hyds_HESS + f1->hyds(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  sucs_HESS   = sucs_HESS + f1->sucs(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  css_HESS    = css_HESS + f1->css(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  rhosoil_HESS= rhosoil_HESS + f1->rhosoil(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
  froot_HESS  = froot_HESS + f1->froot(soil_layer,{A:B},{C:D})*Soil_thickness(soil_layer)/sum(Soil_thickness)
end do

filename2 = path3+"/cable_out_"+year_s+".nc"
f2        = addfile(filename2,"r")

; Albedo_HESS_def = f2->Albedo(0,{A:B},{C:D}) ; Jan

; ================================== Regrid ====================================
; Albedo_LIS_regrid    = regrid_data(Albedo_LIS,interp_method, lat2d, lon2d, \
;               lat, lon, src_mask, dst_mask )
Landcover_LIS_regrid = regrid_data(Landcover_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
Soiltype_LIS_regrid  = regrid_data(Soiltype_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
sand_LIS_regrid      = regrid_data(sand_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
clay_LIS_regrid      = regrid_data(clay_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
silt_LIS_regrid      = regrid_data(silt_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
ssat_LIS_regrid      = regrid_data(ssat_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
bch_LIS_regrid       = regrid_data(bch_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
sfc_LIS_regrid       = regrid_data(sfc_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
swilt_LIS_regrid     = regrid_data(swilt_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
hyds_LIS_regrid      = regrid_data(hyds_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
sucs_LIS_regrid      = regrid_data(sucs_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
css_LIS_regrid       = regrid_data(css_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
rhosoil_LIS_regrid   = regrid_data(rhosoil_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
froot_LIS_regrid     = regrid_data(froot_LIS,interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )
elev_LIS_regrid      = regrid_data(elev_LIS,interp_method, lat2d, lon2d, \
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
pic1 = "./plots/LIS_land_para_output"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"amwg") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")   "amwg") ;

; ; Albedo
; plots = new(3,graphic)
; pres@txString = "Albedo"
; levels = fspan(0,0.3,21)
; res@cnLevels = levels
;
; res@tiMainString = "LIS-CABLE"
; plots(0) = gsn_csm_contour(wks1,Albedo_LIS_regrid*1.,res)
; res@tiMainString = "CABLE(HESS)"
; plots(1) = gsn_csm_contour(wks1,Albedo_HESS*1.,res)
; res@tiMainString = "CABLE(HESS)-def"
; plots(2) = gsn_csm_contour(wks1,Albedo_HESS_def*1.,res)
; gsn_panel(wks1,(/plots/),(/1,3/),pres)
; delete(plots)
; delete(levels)
; delete(res@cnLevels)

; Landcover
plots = new(2,graphic)
pres@txString = "Landcover (mm/y)"
levels = fspan(1,21,21)
res@cnLevels = levels
print(min(Landcover_LIS))
print(max(Landcover_LIS))
res@tiMainString = "LIS-CABLE"
delete(levels)
delete(res@cnLevels)

; ; forest
; tmp_lis = where((Landcover_LIS_regrid .le. 5) .and. (Landcover_LIS_regrid .gt. 0), 1, Landcover_LIS_regrid)
; ; shrub
; tmp_lis = where((tmp_lis .le. 7) .and. (tmp_lis .gt. 5), 2, tmp_lis)
; ; savanna
; tmp_lis = where((tmp_lis .eq. 9) .or. (tmp_lis .eq. 8), 3, tmp_lis)
; ; grassland
; tmp_lis = where(tmp_lis .eq. 10, 4, tmp_lis)
; ; wetland
; tmp_lis = where(tmp_lis .eq. 11, 5, tmp_lis)
; ; cropland
; tmp_lis = where((tmp_lis .eq. 12) .or. (tmp_lis .eq. 14), 6, tmp_lis)
; ; urban
; tmp_lis = where(tmp_lis .eq. 13 , 7, tmp_lis)
; ; ice
; tmp_lis = where(tmp_lis .eq. 15 .or. tmp_lis .eq. 17, 8, tmp_lis)
; ;barren
; tmp_lis = where(tmp_lis .eq. 16 , 9, tmp_lis)
;
; ; forest
; tmp_hess = where((iveg_HESS .le. 4) .and. (iveg_HESS .gt. 0), 1, iveg_HESS)
; ; shrub
; tmp_hess = where(tmp_hess .eq. 5, 2, tmp_hess)
; ; grassland
; tmp_hess = where((tmp_hess .ge. 6) .and. (tmp_hess .le. 8), 4, tmp_hess)
; ; wetland
; tmp_hess = where(tmp_hess .eq. 11, 5, tmp_hess)
; ; cropland
; tmp_hess = where((tmp_hess .eq. 9) .or. (tmp_hess .eq. 10), 6, tmp_hess)
; ; urban
; tmp_hess = where(tmp_hess .eq. 15 , 7, tmp_hess)
; ; ice
; tmp_hess = where(tmp_hess .eq. 17 , 8, tmp_hess)
; ;barren
; tmp_hess = where(tmp_hess .eq. 14 , 9, tmp_hess)
; ;empty
; tmp_hess = where((tmp_hess .eq. 12) .or. (tmp_hess .eq. 13), 10, tmp_hess)

levels = fspan(1,17,17)
res@cnLevels = levels
; plots(0) = gsn_csm_contour(wks1,tmp_lis*1.,res)
plots(0) = gsn_csm_contour(wks1,Landcover_LIS_regrid*1.,res)
res@tiMainString = "CABLE(HESS)"
print(min(iveg_HESS))
print(max(iveg_HESS))
; plots(1) = gsn_csm_contour(wks1,tmp_hess*1.,res)
plots(1) = gsn_csm_contour(wks1,iveg_HESS*1.,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)
delete(levels)
delete(res@cnLevels)

; Soil type
pres@txString = "Soil type"
levels = fspan(1,12,12)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
print(min(Soiltype_LIS))
print(max(Soiltype_LIS))
plots(0) = gsn_csm_contour(wks1,Soiltype_LIS_regrid*1.,res)
res@tiMainString = "CABLE(HESS)"
print(min(isoil_HESS))
print(max(isoil_HESS))
plots(1) = gsn_csm_contour(wks1,isoil_HESS*1.,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

delete(levels)
delete(res@cnLevels)

; bch
pres@txString = "bch"
levels = fspan(0,15.,16)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,bch_LIS_regrid*1.,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,bch_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

delete(levels)
delete(res@cnLevels)

; Sand
pres@txString = "sand"
levels = fspan(0,1.,21)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,sand_LIS_regrid*1.,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,sand_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

; Silt
pres@txString = "silt"
levels = fspan(0,1.,21)
res@cnLevels = levels

res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,silt_LIS_regrid*1.,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,silt_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

; Clay
pres@txString = "clay"
levels = fspan(0,1.,21)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,clay_LIS_regrid*1.,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,clay_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

; ssat
pres@txString = "ssat"
levels = fspan(0,0.5,21)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,ssat_LIS_regrid,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,ssat_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

; sfc
pres@txString = "sfc"
levels = fspan(0,0.5,21)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,sfc_LIS_regrid*1.,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,sfc_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)

; swilt
pres@txString = "swilt"
levels = fspan(0,0.5,21)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,swilt_LIS_regrid*1.,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,swilt_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)
delete(levels)
delete(res@cnLevels)

; hyds
pres@txString = "hyds (m/day)"
levels = fspan(0,2.,21)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,hyds_LIS_regrid*3600*24,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,hyds_HESS/1000.*3600*24,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)
delete(levels)
delete(res@cnLevels)

; sucs
pres@txString = "sucs (m)"
levels = fspan(0,0.5,21)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,abs(sucs_LIS_regrid)*1.,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,sucs_HESS/1000.,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)
delete(levels)
delete(res@cnLevels)

; css
pres@txString = "css"
levels = fspan(700.,1600.,19)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,css_LIS_regrid*1.,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,css_HESS,res)
gsn_panel(wks1,(/plots/),(/1,2/),pres)
delete(levels)
delete(res@cnLevels)
delete(plots)

; rhosoil
plots = new(3,graphic)
filename3 = path4+"/gridinfo_mmy_MD_elev_orig_std_avg-sand_mask.nc"
f3        = addfile (filename3,"r")
rhosoil_grid = (f3->rhosoil_vec(0,{A:B},{C:D})*Soil_thickness(0)\
               +f3->rhosoil_vec(1,{A:B},{C:D})*Soil_thickness(1)\
               +f3->rhosoil_vec(2,{A:B},{C:D})*Soil_thickness(2)\
               +f3->rhosoil_vec(3,{A:B},{C:D})*Soil_thickness(3)\
               +f3->rhosoil_vec(4,{A:B},{C:D})*Soil_thickness(4)\
               +f3->rhosoil_vec(5,{A:B},{C:D})*Soil_thickness(5))/sum(Soil_thickness)

pres@txString = "bulk density (g/cm3)"
levels = fspan(0,3,31)
res@cnLevels = levels
res@tiMainString = "LIS-CABLE"
plots(0) = gsn_csm_contour(wks1,rhosoil_LIS_regrid/1000.,res)
res@tiMainString = "CABLE(HESS)"
plots(1) = gsn_csm_contour(wks1,rhosoil_HESS*(1-ssat_HESS)/1000.,res)
res@tiMainString = "gridinfo"
plots(2) = gsn_csm_contour(wks1,rhosoil_grid/1000.,res)
gsn_panel(wks1,(/plots/),(/1,3/),pres)
; delete(levels)
; delete(res@cnLevels)

; ; froot
; pres@txString = "froot"
; levels = fspan(0,0.2,21)
; res@cnLevels = levels
; res@tiMainString = "LIS-CABLE"
; plots(0) = gsn_csm_contour(wks1,froot_LIS_regrid,res)
; res@tiMainString = "CABLE(HESS)"
; plots(1) = gsn_csm_contour(wks1,froot_HESS,res)
; gsn_panel(wks1,(/plots/),(/1,2/),pres)

delete(plots)
delete(levels)
delete(res@cnLevels)

end
