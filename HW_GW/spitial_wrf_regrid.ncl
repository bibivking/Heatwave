;======================================================================
; ESMF_regrid_5.ncl
;
; Concepts illustrated:
;   - Interpolating from one grid to another using ESMF_regrid
;   - Interpolating data from a WRF grid to a rectilinear grid
;   - Using functions for cleaner code
;======================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;======================================================================
; This procedure plots both the original data and regridded data
; as a panel plot, given the data, variable name, level(s) to
; plot, and the contour levels.
;======================================================================

undef("regrid_data")
function regrid_data(var,Opt,lat2d,lon2d,lat,lon)
local var_regrid
begin
  var@lat2d = lat2d              ; This information will be used by
  var@lon2d = lon2d              ; ESMF_regrid for the source grid
  var_regrid = ESMF_regrid(var,Opt)     ; Do the regridding for TMP
  printVarSummary(var_regrid)
  return(var_regrid)
end

begin

; =================================== Path =====================================
path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/nuwrf_25km/LISWRF_configs/LIS_output"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/outputs"

; ================================= Predefine ==================================
year_s = 2003
year_e = 2008
year_sum = year_e-year_s+1

; study region ??? should larger than SE Aus
A = -43.4377
B = -13.5853
C = 117.552
D = 153.598

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; ================================ READ CABLE ==================================
filename = path1+"/LIS.CABLE."+year_s+"010100.d01.nc"

;---Retrieve either one level, or all levels. Use '-1' for all.
f = addfile(filename,"r")


; ================================== Regrid ====================================
lat2d = f->lat(0,:,:)
lon2d = f->lon(0,:,:)

; study region ??? should larger than SE Aus
A = min(lat2d)
B = max(lat2d)
C = min(lon2d)
D = max(lon2d)

dims  = dimsizes(lat2d)
nlat  = dims(0)
nlon  = dims(1)
lon2d = where(lon2d.lt.0,360+lon2d,lon2d)
;---Create the destination lat/lon grid
lat = fspan(A,B,nlat)
lon = fspan(C,D,nlon)

Opt                = True
Opt@SrcTitle       = "WRF grid"   ; optional
Opt@ForceOverwrite = True

Opt@SrcFileName    = "WRF_SCRIP.nc"      ; Name of source and ???
Opt@DstFileName    = "Rectilinear.nc"    ; destination files  ???


; Opt@SrcGridCornerLat = sfile->lat_vertices    ; corners are necessary
; Opt@SrcGridCornerLon = sfile->lon_vertices    ; for "conserve" method
Opt@SrcGridMask      = where(.not.ismissing(f->Evap_tavg(0,:,:)),1,0)

; Opt@DstGridType    = "rectilinear"
; Opt@DstGridLat     = lat
; Opt@DstGridLon     = lon

Opt@DstGridType      = "1deg"              ; Destination grid
; Opt@DstTitle         = "1-degree Resolution"
; Opt@DstLLCorner      = (/-43,  117 /)
; Opt@DstURCorner      = (/-14,  154 /)

Opt@InterpMethod   = "neareststod"; "conserve" "patch";"bilinear"
Opt@SrcRegional    = True
Opt@DstRegional    = True

;Opt@DstGridMask    = where(.not.ismissing(thetao),1,0)




; if first_time /= True then
;   Opt@SkipSrcGrid   = True
;   Opt@SkipDstGrid   = True
;   Opt@SkipWgtGen    = True
; end if


var_regrid = regrid_data(f->Evap_tavg(0,:,:),Opt,lat2d,lon2d,lat,lon)

; _________________________ Check Value ____________________________
pic = "check_regrid"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
res            = True
res@cnFillMode = "RasterFill"            ; Raster Mode
res@cnFillOn   = True
res@cnLinesOn  = False

plot = gsn_csm_contour(wks,f->Evap_tavg(0,:,:),res)
plot = gsn_csm_contour(wks,var_regrid,res)
end
