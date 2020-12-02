;======================================================================
; ESMF_regrid_32.ncl
;======================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

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

    ;Opt@DstGridType       = "1.0deg"

    ; if first_time /= True then
    ;   Opt@SkipSrcGrid   = True
    ;   Opt@SkipDstGrid   = True
    ;   Opt@SkipWgtGen    = True
    ; end if

    var@lat2d = lat2d              ; This information will be used by
    var@lon2d = lon2d              ; ESMF_regrid for the source grid
    var_regrid = ESMF_regrid(var,Opt)     ; Do the regridding for TMP
    printVarSummary(var_regrid)
    return(var_regrid)
end

begin

; =================================== Path =====================================
path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/nuwrf_25km/LISWRF_configs/LIS_output"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/outputs/GSWP"

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
f = addfile(filename,"r")

filename1 = path2+"/cable_out_1979.nc"
f2 = addfile(filename1,"r")
; ================================== Regrid ====================================
interp_method = "neareststod"  ;"neareststod"; "conserve" "patch";"bilinear"

minlon  = 117 ;min(lon2d)-1 ; 116.5516
minlat  = -44 ;min(lat2d)-1 ; -44.43771
maxlon  = 154 ;max(lon2d)+1 ; 154.5981
maxlat  = -13 ;max(lat2d)+1 ; -12.5853

src_mask = where(.not.ismissing(f->Evap_tavg(0,:,:)),1,0)
dst_mask = where(.not.ismissing(f2->Qs(0,{minlat:maxlat},{minlon:maxlon})),1,0)

lat2d    = f->lat(0,:,:)
lon2d    = f->lon(0,:,:)
lon2d    = where(lon2d.lt.0,360+lon2d,lon2d)
print(dimsizes(dst_mask))
lat      = fspan(minlat,maxlat,dimsizes(dst_mask(:,0)))
lon      = fspan(minlon,maxlon,dimsizes(dst_mask(0,:)))

var_regrid = regrid_data(f->Evap_tavg(0,:,:),interp_method, lat2d, lon2d, \
              lat, lon, src_mask, dst_mask )

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
plot = gsn_csm_contour(wks,src_mask,res)
plot = gsn_csm_contour(wks,dst_mask,res)

end
