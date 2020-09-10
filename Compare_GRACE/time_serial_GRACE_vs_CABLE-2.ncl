;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

undef("regrid_data")
function regrid_data(var,interp_method,lat2d,lon2d,lat,lon,src_mask,dst_mask)
local var_regrid
begin

    Opt                   = True

    Opt@SrcFileName       = "./nc_files/WRF_SCRIP.nc"      ; Name of source and ???
    Opt@DstFileName       = "./nc_files/1degree.nc"    ; destination files  ???

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

; ================================= Start ======================================
begin

case_name = (/"spinup2"/);(/"spinup1","spinup2"/)
      ;(/ "Princeton_ctl_watmove_new_veg" /)
      ; (/ "GSWP" /)
path = "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3/outputs/GSWP/"
      ; "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"
      ; "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3/outputs/GSWP/"
      ; "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"
      ; "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3/outputs/"

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

; =============================== Pre-define ===================================
year_s   = 1979
year_e   = 2010
GRACE_start_year = 2003
GRACE_end_year   = 2009
year_sum = year_e-year_s+1
total    = year_sum*12
case_sum = dimsizes(case_name)


; ================================ Set Regrid ==================================
use_month     = "May" ; "Jul"
time_scale    = "monthly" ; "annual"; "monthly"
grace_version = "Gap-fill" ;"Gap-fill"; "Original"
output        = "daily"
cable_version = "CABLE-HESS" ;"LIS-CABLE" ; "CABLE-HESS"

if cable_version .eq. "LIS-CABLE" then
    scale = "AU-domain"
    interp_method = "neareststod"  ;"neareststod"; "conserve" "patch";"bilinear"

    filename  = path+case_name(0)+"/LIS_output/LIS.CABLE."+year_s+"090100.d01.nc"
    f = addfile (filename,"r")

    lat_size = dimsizes(f->lat(0,:,0))
    lon_size = dimsizes(f->lat(0,0,:))
    default_value = -9999.
    lat2d    = f->lat(0,:,:)
    lon2d    = f->lon(0,:,:)
    lon2d    = where(lon2d.lt.0,360+lon2d,lon2d)

    minlon  = ceil(min(lon2d))+0.5 ; 116.5516
    minlat  = ceil(min(lat2d))+0.5 ; -44.43771
    maxlon  = floor(max(lon2d))-0.5 ; 154.5981
    maxlat  = floor(max(lat2d))-0.5 ; -12.5853

    ; study region
    A = minlat
    B = maxlat
    C = minlon
    D = maxlon
    delete(f)
else if cable_version .eq. "CABLE-HESS" then
    scale = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"

    filename = path+case_name(0)+"/cable_out_"+year_s+".nc"
    f = addfile (filename,"r")

    lat_size  = dimsizes(f->y)
    lon_size  = dimsizes(f->x)
    CABLE_lat = f->y
    CABLE_lon = f->x
    default_value = -1.e+33

    if scale .eq. "Global" then
      A = -90.
      B = 90.
      C = 0.
      D = 360.
    else if scale .eq. "SE-AU" then
      A = -40.
      B = -27.
      C = 140.
      D = 155.
    else if scale .eq. "AU" then
      A = -44.
      B = -10.
      C = 112.
      D = 154.
    end if
    end if
    end if
    delete(f)
end if
end if

; ================================= Check ======================================
pic = "check"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
res = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;???
res@tmBorderThicknessF = 3.0
res@cnLinesOn          = False                       ; ???????????
res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ?????

; =========================== Read GRACE and fill gap ==========================
GRACE      = new((/total,360,720/),float)
uncertainty= new((/total,360,720/),float)
GRACE      = 0.
uncertainty= 0.

if GRACE_start_year .gt. year_s then
    compensate_month = (GRACE_start_year-year_s) * 12
    ; months before 2003
else
    compensate_month = 0
end if

if grace_version .eq. "Original" then
    filename   = "/g/data/w35/mm3972/data/GRACE/GRACE_JPL_RL06/GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc"
    f          = addfile (filename,"r")

    GRACE(compensate_month+0:compensate_month+4,:,:)       = dble2flt(f->lwe_thickness(7:11,:,:))*10.
    GRACE(compensate_month+6:compensate_month+total-1,:,:) = dble2flt(f->lwe_thickness(12:total+5,:,:))*10.

    uncertainty(compensate_month+0:compensate_month+4,:,:) = dble2flt(f->uncertainty(7:11,:,:))*10.
    uncertainty(compensate_month+6:compensate_month+total-1,:,:)= dble2flt(f->uncertainty(12:total+5,:,:))*10.

    printVarSummary(GRACE)
    printVarSummary(uncertainty)

    May_ave = new((/360,720/),float)
    Jun_ave = new((/360,720/),float)
    Jul_ave = new((/360,720/),float)
    May_ave = 0.
    Jun_ave = 0.
    Jul_ave = 0.

    GRACE@_FillValue       = dble2flt(f->lwe_thickness@_FillValue)
    uncertainty@_FillValue = dble2flt(f->lwe_thickness@_FillValue)

    printVarSummary(GRACE)
    printVarSummary(uncertainty)

    plot = gsn_csm_contour(wks,GRACE(compensate_month,:,:),res)

    do year = GRACE_start_year+1,GRACE_end_year
        print("add year :" + year)
        May_ave = May_ave + GRACE(compensate_month+(year-GRACE_start_year)*12+4,:,:)
        Jun_ave = Jun_ave + GRACE(compensate_month+(year-GRACE_start_year)*12+5,:,:)
        Jul_ave = Jul_ave + GRACE(compensate_month+(year-GRACE_start_year)*12+6,:,:)
    end do

    May_ave = May_ave/(GRACE_end_year-GRACE_start_year)
    Jun_ave = Jun_ave/(GRACE_end_year-GRACE_start_year)
    Jul_ave = Jul_ave/(GRACE_end_year-GRACE_start_year)
    print(May_ave)
    plot = gsn_csm_contour(wks,May_ave,res)
    plot = gsn_csm_contour(wks,Jun_ave,res)
    plot = gsn_csm_contour(wks,Jul_ave,res)

    if use_month .eq. "Jul" then
       GRACE(compensate_month+5,:,:) = GRACE(compensate_month+6,:,:) * (Jun_ave/Jul_ave)
    else if use_month .eq. "May" then
       GRACE(compensate_month+5,:,:) = GRACE(compensate_month+4,:,:) * (Jun_ave/May_ave)
    end if
    end if

    uncertainty(compensate_month+5,:,:) = (uncertainty(compensate_month+4,:,:) + uncertainty(compensate_month+6,:,:))/2.

    delete(May_ave)
    delete(Jun_ave)
    delete(Jul_ave)

    plot = gsn_csm_contour(wks,GRACE(compensate_month+5,:,:),res)
    plot = gsn_csm_contour(wks,uncertainty(compensate_month+5,:,:),res)

else if grace_version .eq. "Gap-fill" then

    if GRACE_start_year .lt. year_s then
      GRACE_start_year = year_s
    end if

    if GRACE_end_year .gt. year_e then
      GRACE_end_year = year_e
    end if

    do year = GRACE_start_year,GRACE_end_year

        filename   = "/g/data/w35/mm3972/data/GRACE/GRACE_JPL_RL06/fill_gap_filter_singular_Sanaa/GRACERel06_lwe_"+ year+".nc"
        f          = addfile (filename,"r")

        a = (year-year_s)*12
        b = (year-year_s)*12+11

        GRACE(a:b,:,:)       = lonFlip(f->lwe)*10.
        uncertainty(a:b,:,:) = lonFlip(f->uncertainty)*10.

    end do

    printVarSummary(GRACE)
    printVarSummary(uncertainty)

    GRACE@_FillValue       = -1.e+33
    uncertainty@_FillValue = -1.e+33

    replace_ieeenan (GRACE, GRACE@_FillValue, 0)
    replace_ieeenan (uncertainty, uncertainty@_FillValue, 0)

    printVarSummary(GRACE)
    printVarSummary(uncertainty)

    plot = gsn_csm_contour(wks,GRACE(0,:,:),res)
    plot = gsn_csm_contour(wks,uncertainty(0,:,:),res)

end if
end if


; ====== post process GRACE =======
filename = "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3/outputs/GSWP/spinup2/cable_out_1999.nc"
f = addfile (filename,"r")
landmask = where(ismissing(f->Qle(0,:,:)), 0, 1)

GRACE_lat = f->y
GRACE_lon = f->x

GRACE!0 = "time"
GRACE!1 = "lat"
GRACE!2 = "lon"
GRACE&lat = GRACE_lat
GRACE&lon = GRACE_lon

uncertainty!0 = "time"
uncertainty!1 = "lat"
uncertainty!2 = "lon"
uncertainty&lat = GRACE_lat
uncertainty&lon = GRACE_lon

printVarSummary(landmask)
do time = 0, total-1
  GRACE(time,:,:)       = where( landmask .eq. 1, GRACE(time,:,:), GRACE@_FillValue)
  uncertainty(time,:,:) = where( landmask .eq. 1, uncertainty(time,:,:), uncertainty@_FillValue)
end do
printVarSummary(GRACE)
printVarSummary(uncertainty)

delete(f)
delete(filename)
delete(landmask)

;================================ Read CABLE ===================================
TWS           = new((/case_sum,total,lat_size, lon_size/),float)
GWMoist       = new((/case_sum,total,lat_size, lon_size/),float)
SoilMoist     = new((/case_sum,total,lat_size, lon_size/),float)
CanopInt      = new((/case_sum,total,lat_size, lon_size/),float)
SWE           = new((/case_sum,total,lat_size, lon_size/),float)
Rain          = new((/case_sum,total,lat_size, lon_size/),float)

GWMoist       = 0.0
SoilMoist     = 0.0
CanopInt      = 0.0
SWE           = 0.0
TWS           = 0.0
Rain          = 0.0

if cable_version .eq. "CABLE-HESS" then
    do case_num = 0,case_sum-1
       print(case_name(case_num))
       counter = 0
       do year = year_s,year_e
          print("year = " +year)
          filename1 = path+case_name(case_num)+"/cable_out_"+year+".nc"
          f1 = addfile (filename1,"r")
          if output .eq. "monthly" then
            do month = 0,11
               TWS(case_num,counter,:,:)       = f1->GWMoist(month,:,:)*f1->GWdz*1000. + f1->CanopInt(month,:,:) + f1->SWE(month,:,:)
               GWMoist(case_num,counter,:,:)   = f1->GWMoist(month,:,:)*f1->GWdz*1000.
               CanopInt(case_num,counter,:,:)  = f1->CanopInt(month,:,:)
               SWE(case_num,counter,:,:)       = f1->SWE(month,:,:)
               Rain(case_num,counter,:,:)      = f1->Rainf(month,:,:)
               do layer = 0,5
                  SoilMoist(case_num,counter,:,:) = SoilMoist(case_num,counter,:,:) + f1->SoilMoist(month,layer,:,:)*Soil_thickness(layer)*1000.
                  TWS(case_num,counter,:,:)       = TWS(case_num,counter,:,:) + f1->SoilMoist(month,layer,:,:)*Soil_thickness(layer)*1000.
               end do
               counter = counter + 1
            end do ;month
          else if output .eq. "daily" then
            if ( mod(year,4) .eq. 0) then
               dom = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
            else
               dom = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
            end if
            do month = 0,11
               TWS(case_num,counter,:,:) = dim_avg_n_Wrap(f1->GWMoist(dom(month):dom(month+1)-1,:,:),0) \
                                          *f1->GWdz*1000. + dim_avg_n_Wrap(f1->CanopInt(dom(month):dom(month+1)-1,:,:),0)\
                                          + dim_avg_n_Wrap(f1->SWE(dom(month):dom(month+1)-1,:,:), 0)
               GWMoist(case_num,counter,:,:)   = dim_avg_n_Wrap(f1->GWMoist(dom(month):dom(month+1)-1,:,:),0)*f1->GWdz*1000.
               CanopInt(case_num,counter,:,:)  = dim_avg_n_Wrap(f1->CanopInt(dom(month):dom(month+1)-1,:,:),0)
               SWE(case_num,counter,:,:)       = dim_avg_n_Wrap(f1->SWE(dom(month):dom(month+1)-1,:,:), 0)
               Rain(case_num,counter,:,:)      = dim_avg_n_Wrap(f1->Rainf(dom(month):dom(month+1)-1,:,:), 0)
               do layer = 0,5
                  SoilMoist(case_num,counter,:,:) = SoilMoist(case_num,counter,:,:) + dim_avg_n_Wrap( \
                                             f1->SoilMoist(dom(month):dom(month+1)-1,layer,:,:),0)\
                                             *Soil_thickness(layer)*1000.
                  TWS(case_num,counter,:,:) = TWS(case_num,counter,:,:) + dim_avg_n_Wrap( \
                                             f1->SoilMoist(dom(month):dom(month+1)-1,layer,:,:),0)\
                                             *Soil_thickness(layer)*1000.
               end do
               counter = counter + 1
            end do ;month
          end if ; daily
          end if
       end do ;year
    end do ;case

    TWS!0 = "case"
    TWS!1 = "time"
    TWS!2 = "lat"
    TWS!3 = "lon"
    TWS&lat = CABLE_lat
    TWS&lon = CABLE_lon

    GWMoist!0 = "case"
    GWMoist!1 = "time"
    GWMoist!2 = "lat"
    GWMoist!3 = "lon"
    GWMoist&lat = CABLE_lat
    GWMoist&lon = CABLE_lon

    CanopInt!0 = "case"
    CanopInt!1 = "time"
    CanopInt!2 = "lat"
    CanopInt!3 = "lon"
    CanopInt&lat = CABLE_lat
    CanopInt&lon = CABLE_lon

    SWE!0 = "case"
    SWE!1 = "time"
    SWE!2 = "lat"
    SWE!3 = "lon"
    SWE&lat = CABLE_lat
    SWE&lon = CABLE_lon

    SoilMoist!0 = "case"
    SoilMoist!1 = "time"
    SoilMoist!2 = "lat"
    SoilMoist!3 = "lon"
    SoilMoist&lat = CABLE_lat
    SoilMoist&lon = CABLE_lon

    Rain!0 = "case"
    Rain!1 = "time"
    Rain!2 = "lat"
    Rain!3 = "lon"
    Rain&lat = CABLE_lat
    Rain&lon = CABLE_lon
else if cable_version .eq. "LIS-CABLE" then

    GWdz = 25. ; aquifer thickness in LIS is fixed to 25m

    do case_num = 0,case_sum-1
       print(case_name(case_num))
       counter = 0
       do year = year_s,year_e
          print("year = " +year)
          do month = 1, 12
              if month .lt. 10 then
                filename1  = path+case_name(case_num)+"/LIS_output/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
              else
                filename1  = path+case_name(case_num)+"/LIS_output/LIS.CABLE."+year+month+"0100.d01.nc"
              end if
              f1  = addfile (filename1,"r")
              TWS(case_num,counter,:,:) = dim_avg_n_Wrap(f1->GWwb_tavg(:,:,:),0)*GWdz*1000. \
                                        + dim_avg_n_Wrap(f1->CanopInt_inst(:,:,:),0)\
                                        + dim_avg_n_Wrap(f1->SWE_inst(:,:,:), 0)
              GWMoist(case_num,counter,:,:)   = dim_avg_n_Wrap(f1->GWwb_tavg(:,:,:),0)*GWdz*1000.
              CanopInt(case_num,counter,:,:)  = dim_avg_n_Wrap(f1->CanopInt_inst(:,:,:),0)
              SWE(case_num,counter,:,:)       = dim_avg_n_Wrap(f1->SWE_inst(:,:,:), 0)
              Rain(case_num,counter,:,:)      = dim_avg_n_Wrap(f1->Rainf_tavg(:,:,:), 0)

              do layer = 0,5
                  SoilMoist(case_num,counter,:,:) = SoilMoist(case_num,counter,:,:) +  \
                                             dim_avg_n_Wrap(f1->SoilMoist_tavg(:,layer,:,:),0)\
                                             *Soil_thickness(layer)*1000.
                  TWS(case_num,counter,:,:) = TWS(case_num,counter,:,:) +  \
                                              dim_avg_n_Wrap(f1->SoilMoist_tavg(:,layer,:,:),0)\
                                              *Soil_thickness(layer)*1000.
              end do
              counter = counter + 1
          delete(f1)
          end do ; month
       end do ; year
    end do ; case
end if ; cable_version
end if

TWS@_FillValue       = default_value
GWMoist@_FillValue   = default_value
CanopInt@_FillValue  = default_value
SWE@_FillValue       = default_value
SoilMoist@_FillValue = default_value
Rain@_FillValue      = default_value

plot = gsn_csm_contour(wks,GRACE(compensate_month,:,:),res)
plot = gsn_csm_contour(wks,uncertainty(compensate_month,:,:),res)
plot = gsn_csm_contour(wks,TWS(0,5,:,:),res)

; Rain_LIS_regrid   = regrid_data(Rain_LIS,interp_method, lat2d, lon2d, \
;                         lat, lon, src_mask, dst_mask )


; =========================== Calculate 2003-2009 avg ==========================
TWS_ave       = new((/case_sum, lat_size, lon_size/),float)
GWMoist_ave   = new((/case_sum, lat_size, lon_size/),float)
SoilMoist_ave = new((/case_sum, lat_size, lon_size/),float)
CanopInt_ave  = new((/case_sum, lat_size, lon_size/),float)
SWE_ave       = new((/case_sum, lat_size, lon_size/),float)
TWS_ave       = 0.0
GWMoist_ave   = 0.0
SoilMoist_ave = 0.0
CanopInt_ave  = 0.0
SWE_ave       = 0.0

; Calculate 2004-2009 TWS average
if year_e .gt. 2009 then
  adjust  = (year_e - GRACE_end_year)*12
else
  adjust  = 0
  GRACE_end_year = year_e
end if
TWS_ave        =  dim_avg_n_Wrap(TWS(:,compensate_month+12:total-1 - adjust,:,:),1)
GWMoist_ave    =  dim_avg_n_Wrap(GWMoist(:,compensate_month+12:total-1 - adjust,:,:),1)
CanopInt_ave   =  dim_avg_n_Wrap(CanopInt(:,compensate_month+12:total-1 - adjust,:,:),1)
SWE_ave        =  dim_avg_n_Wrap(SWE(:,compensate_month+12:total-1 - adjust,:,:),1)
SoilMoist_ave  =  dim_avg_n_Wrap(SoilMoist(:,compensate_month+12:total-1 - adjust,:,:),1)

; Calculate TWS anormaly
do time = 0,total -1
    TWS(:,time,:,:)      = TWS(:,time,:,:) - TWS_ave
    GWMoist(:,time,:,:)  = GWMoist(:,time,:,:) - GWMoist_ave
    CanopInt(:,time,:,:) = CanopInt(:,time,:,:) - CanopInt_ave
    SWE(:,time,:,:)      = SWE(:,time,:,:) - SWE_ave
    SoilMoist(:,time,:,:)= SoilMoist(:,time,:,:) - SoilMoist_ave
end do

delete(TWS_ave)
delete(GWMoist_ave)
delete(CanopInt_ave)
delete(SWE_ave)
delete(SoilMoist_ave)

; =============================  Regional AVE  =================================
print("global average")
GRACE_global_ave       = new((/3,total/),float)
uncertainty_global_ave = new(total,float)
TWS_global_ave         = new((/case_sum,total/),float)
GWMoist_global_ave     = new((/case_sum,total/),float)
CanopInt_global_ave    = new((/case_sum,total/),float)
SWE_global_ave         = new((/case_sum,total/),float)
SoilMoist_global_ave   = new((/case_sum,total/),float)
Rain_global_ave        = new((/case_sum,total/),float)

GRACE_global_ave       = 0.0
uncertainty_global_ave = 0.0
TWS_global_ave         = 0.0
GWMoist_global_ave     = 0.0
CanopInt_global_ave    = 0.0
SWE_global_ave         = 0.0
SoilMoist_global_ave   = 0.0
Rain_global_ave        = 0.0

do time = 0,total-1
   GRACE_global_ave(0,time) = avg(GRACE(time,{A:B},{C:D}))
   uncertainty_global_ave(time) = avg(uncertainty(time,{A:B},{C:D}))
   do case_num = 0,case_sum-1
      if cable_version .eq. "CABLE-HESS" then
          TWS_global_ave(case_num,time)       = avg(TWS(case_num,time,{A:B},{C:D}))
          GWMoist_global_ave(case_num,time)   = avg(GWMoist(case_num,time,{A:B},{C:D}))
          CanopInt_global_ave(case_num,time)  = avg(CanopInt(case_num,time,{A:B},{C:D}))
          SWE_global_ave(case_num,time)       = avg(SWE(case_num,time,{A:B},{C:D}))
          SoilMoist_global_ave(case_num,time) = avg(SoilMoist(case_num,time,{A:B},{C:D}))
          Rain_global_ave(case_num,time)      = avg(Rain(case_num,time,{A:B},{C:D}))
      else if cable_version .eq. "LIS-CABLE" then
          TWS_global_ave(case_num,time)       = avg(TWS(case_num,time,:,:))
          GWMoist_global_ave(case_num,time)   = avg(GWMoist(case_num,time,:,:))
          CanopInt_global_ave(case_num,time)  = avg(CanopInt(case_num,time,:,:))
          SWE_global_ave(case_num,time)       = avg(SWE(case_num,time,:,:))
          SoilMoist_global_ave(case_num,time) = avg(SoilMoist(case_num,time,:,:))
          Rain_global_ave(case_num,time)      = avg(Rain(case_num,time,:,:))
      end if
      end if
   end do
end do

GRACE_global_ave(1,:) = GRACE_global_ave(0,:) - uncertainty_global_ave
GRACE_global_ave(2,:) = GRACE_global_ave(0,:) + uncertainty_global_ave

GRACE_TWS              = new(case_sum,float)
Rain_TWS               = new(case_sum,float)
GRACE_TWS              = 0.0
Rain_TWS               = 0.0

Rain_GRACE = escorc(GRACE_global_ave(0,:), Rain_global_ave(0,:))
print("Correlation between Rain and GRACE : "+Rain_GRACE)

do case_num = 0,case_sum-1
  GRACE_TWS(case_num) = escorc(GRACE_global_ave(0,:), TWS_global_ave(case_num,:))
  print("Correlation between GRACE and TWS in "+case_name(case_num)+": "+GRACE_TWS(case_num))

  Rain_TWS(case_num) = escorc(Rain_global_ave(0,:), TWS_global_ave(case_num,:))
  print("Correlation between Rain and TWS in "+case_name(case_num)+": "+Rain_TWS(case_num))
end do

GRACE_global_ave_year       = new((/3,year_sum/),float)
uncertainty_global_ave_year = new(year_sum,float)
TWS_global_ave_year         = new((/case_sum,year_sum/),float)
GRACE_global_ave_year       = 0.0
uncertainty_global_ave_year = 0.0
TWS_global_ave_year         = 0.0

do year = year_s,year_e
  GRACE_global_ave_year(:,year-year_s) = dim_avg_n_Wrap(GRACE_global_ave(:,((year-year_s)*12):((year-year_s)*12+11)),1)
  TWS_global_ave_year(:,year-year_s)   = dim_avg_n_Wrap(TWS_global_ave(:,((year-year_s)*12):((year-year_s)*12+11)),1)
  uncertainty_global_ave_year(year-year_s) = dim_avg_n_Wrap(uncertainty_global_ave(((year-year_s)*12):((year-year_s)*12+11)),0)
end do

delete(GRACE)
delete(TWS)
delete(uncertainty)

;GRACE & TWS is anormaly here
delete(wks)
delete(res)
delete(pic)

;_____________________________ PLOT ___________________________________
pic = "./plots/time_serial_GRACE_JPL06-"+grace_version+"_vs_"+cable_version+"-"+case_name(0)+"_"+scale+"_"+time_scale

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

pres                   = True
res                    = True              ; plot mods desired
res@gsnDraw            = False             ; don't draw yet
res@gsnFrame           = False             ; don't advance frame yet


res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
res@vpWidthF 	   = 0.7
res@txFontHeightF	 = 0.015            ; change title font heights

res@trYMaxF = 150.
res@trYMinF = -150.

LineColors       = (/"darkgreen","limegreen","red","pink","orange"/) ; colors chosen
LineThicknesses	 = (/2.,2.,2.,2.,2./)      ; line thicknesses
DashPatterns	   = (/0.,0.,0.,0.,0./)      ; make all lines solid
Labels           = (/"TWS", "GWMoist","SoilMoist", "CanopInt", "SWE" /)

if time_scale .eq. "monthly" then

    x = ispan(1,total,1)
    res@tmXBMode  = "Explicit"
    res@tmXBValues= ispan(1,total+1,12*5)
    res@tmXBLabels= tostring(ispan(year_s,year_e+1,5))

    res@xyLineThicknessF       = 2.
    res@xyLineColors           = LineColors
    res@xyLineThicknesses      = LineThicknesses
    res@xyDashPatterns	       = DashPatterns
    res@xyExplicitLegendLabels = Labels

    do case_num = 0,case_sum - 1
        res@txString           = grace_version+" vs "+cable_version+" "+case_name(case_num)

        res@pmLegendDisplayMode= "Never"
        res@xyLineColor        = "Blue"
        res@xyMonoLineColor    = True             ; want colored lines
        top_plot = gsn_csm_xy (wks,x,GRACE_global_ave(0,:),res)       ; create line plot
        delete(res@xyLineColor)

        ; Create a rain plot
        res@xyLineColor        = "gray"
        bot_plot  = gsn_csm_xy (wks,x,Rain_global_ave(0,:),res)  ; Create filled XY plot.

      ; Create a plot with the area between both curves filled in blue.
      ;  res@gsnXYFillColors = "LightBlue"
      ;  res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
      ;  bot_plot  = gsn_csm_xy (wks,x,GRACE_global_ave(1:2,:),res)  ; Create filled XY plot.
      ;  delete(res@gsnXYFillColors)

      ; Create a plot with the area between both curves filled in pink.


        res@pmLegendDisplayMode    = "Always"
        res@xyMonoLineColor        = False             ; want colored lines
        tmp = new((/5,total/),float)

        tmp(0,:) = TWS_global_ave(case_num,:)
        tmp(1,:) = GWMoist_global_ave(case_num,:)
        tmp(2,:) = SoilMoist_global_ave(case_num,:)
        tmp(3,:) = CanopInt_global_ave(case_num,:)
        tmp(4,:) = SWE_global_ave(case_num,:)

        mid_plot  = gsn_csm_xy (wks,x,tmp,res)  ; Create another filled XY plot.
        delete(tmp)

        overlay(mid_plot,top_plot)
        ; overlay(bot_plot,mid_plot)
        ; overlay(bot_plot,top_plot)

        maximize_output(wks,pres)

        delete(top_plot)
        delete(mid_plot)
        delete(bot_plot)

    end do; case_num
else if time_scale .eq. "annual" then

  x = ispan(1,year_sum,1)
  res@tmXBValues=ispan(1,year_sum,1)
  xlabals = (/"2003","2004","2005","2006","2007","2008","2009","2010"/)
  res@tmXBLabels=xlabals(0:year_sum-1)

  res@xyMonoLineColor    = False             ; want colored lines
  res@xyLineColor        = "Blue"
  res@xyLineThicknessF   = 2.
  top_plot = gsn_csm_xy (wks,x,GRACE_global_ave_year(0,:),res)       ; create line plot

; Create a plot with the area between both curves filled in blue.
 delete(res@xyLineColor)
 res@gsnXYFillColors    = "LightBlue"
 res@xyLineColor        = -1                           ; We don't want the line, so make it transparent.
 bot_plot  = gsn_csm_xy (wks,x,GRACE_global_ave_year(1:2,:),res)  ; Create filled XY plot.

; Create a plot with the area between both curves filled in pink.
  delete(res@gsnXYFillColors)
  res@xyLineColors       = LineColors
  res@xyLineThicknesses	 = LineThicknesses
  res@xyDashPatterns	   = DashPatterns
  mid_plot  = gsn_csm_xy (wks,x,TWS_global_ave_year,res)  ; Create another filled XY plot.

   overlay(mid_plot,top_plot)
;  overlay(bot_plot,mid_plot)
;  overlay(bot_plot,top_plot)
;  maximize_output(wks,pres)

end if
end if

; ;*****************************************************
; ; Manually create and attach legend
; ;*****************************************************
; dum1 = new(case_sum,graphic)
; dum2 = new(case_sum,graphic)
;
; res_text                = True                  ; text mods desired
; res_text@txFontHeightF  = 0.013                ; change text size
; res_text@txJust         = "CenterLeft"          ; text justification
; res_lines               = True                  ; polyline mods desired
;
; do case_num = 0,case_sum-1
;   if case_num .lt. case_sum then
;      res_lines@gsLineDashPattern = DashPatterns(case_num)
;      res_lines@gsLineThicknessF  = LineThicknesses(case_num)
;      res_lines@gsLineColor       = LineColors(case_num)
;   else
;      res_lines@gsLineDashPattern = 0.
;      res_lines@gsLineThicknessF  = 2.
;      res_lines@gsLineColor       = "blue"                                 ; change to blue
;   end if
;
;   if case_num .eq. 7 then
;      point_xl = 1
;      point_xr = 5
;      point_xt = 5
;      point_y  = 65
;   else if case_num .ge. 0 .and. case_num .le. 2 then
;      point_xl = (case_num+1)*24+1
;      point_xr = (case_num+1)*24+6
;      point_xt = (case_num+1)*24+6
;      point_y = 65
;   else
;      point_xl = (case_num-3)*24+1
;      point_xr = (case_num-3)*24+6
;      point_xt = (case_num-3)*24+6
;      point_y = 58
;   end if
;   end if
;   xx = (/point_xl,point_xr/)
;   yy = (/point_y,point_y/)
;   dum1(case_num) = gsn_add_polyline(wks,mid_plot,xx,yy,res_lines)              ; add polyline
;   if case_num .lt. case_sum then
;      dum2(case_num) = gsn_add_text(wks,mid_plot,case_name(case_num),point_xt,point_y,res_text); add text
;   else
;      dum2(case_num) = gsn_add_text(wks,mid_plot,"GRACE",point_xt,point_y,res_text); add text
;   end if
; end do
; maximize_output(wks,pres)

end
