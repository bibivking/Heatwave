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

case_name = (/"FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
            ; "CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr", \
            ; "CTL-30x3+30yr-DeepRoot/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
path      = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
pic_message = "FD+DR"

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

; =============================== Pre-define ===================================
year_s   = 2000
year_e   = 2019
GRACE_start_year = 2002
GRACE_end_year   = year_e ; 2019
year_sum = year_e-year_s+1
total    = year_sum*12
case_sum = dimsizes(case_name)

; ================================ Set Regrid ==================================
PFT           = "tree" ; "tree" ;"all"
use_month     = "May" ; "Jul"
time_scale    = "monthly" ; "annual"; "monthly"
grace_version = "Original-fill-nan"
                ;"Original-fill-nan": using Nan fill missing month
                ;"Gap-fill" : Sanaa's gap fill
                ;"Original" : if this option, should alwasy set GRACE_start_year = 2003
                ;             GRACE_end_year <= 2010, since too many gaps after 2010
output        = "daily" ;"monthly"
cable_version = "CABLE-HESS" ;"LIS-CABLE" ; "CABLE-HESS"

; set PFT
if PFT .eq. "all" then
  iveg_s = 0
  iveg_e = 20
else if PFT .eq. "tree" then
  iveg_s = 1
  iveg_e = 4
end if
end if

; ================================== Set Dim ===================================
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

    filename = path+case_name(0)+"/cable_out_"+year_s+"_SE_Aus.nc" ; MMY
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
      B = -28.
      C = 140.
      D = 154.
    else if scale .eq. "AU" then
      A = -44.
      B = -10.
      C = 112.
      D = 154.
    end if
    end if
    end if

    ; PFT mask
    PFT_mask  = where((f->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
    PFT_mask_multi = new((/total,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
    do z = 0, total-1
      PFT_mask_multi(z,:,:) = PFT_mask
    end do
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
GRACE      = new((/total,360,720/),float, -99999.)
uncertainty= new((/total,360,720/),float, -99999.)
GRACE      = -99999.
uncertainty= -99999.

; compensate_month : months before GRACE start year
if GRACE_start_year .gt. year_s then
    compensate_month = (GRACE_start_year-year_s) * 12
else
    compensate_month = 0
end if

if grace_version .eq. "Original-fill-nan" then

  filename   = "/g/data/w35/Shared_data/Observations/GRACE/GRACE_JPLRL06M_MASCON/GRCTellus.JPL.2002_2020_gap-fill_MMY.nc"
  f          = addfile (filename,"r")

  tot_grace       = (GRACE_end_year-GRACE_start_year+1)*12 ; how many year will be read from GRACE
  GRACE_ts_str_yr = (GRACE_start_year - 2002) * 12
  GRACE_ts_tot    = (GRACE_end_year - GRACE_start_year + 1) * 12 -1
  GRACE(compensate_month:compensate_month+tot_grace-1,:,:)      = \
            f->lwe_thickness(GRACE_ts_str_yr:GRACE_ts_str_yr+GRACE_ts_tot, :, :)*10.
  uncertainty(compensate_month:compensate_month+tot_grace-1,:,:)= \
            f->uncertainty(GRACE_ts_str_yr:GRACE_ts_str_yr+GRACE_ts_tot, :, :)*10.
  replace_ieeenan (GRACE, GRACE@_FillValue, 0)
  replace_ieeenan (uncertainty, uncertainty@_FillValue, 0)

  printVarSummary(GRACE)
  printVarSummary(uncertainty)

else if grace_version .eq. "Original" then
    filename   = "/g/data/w35/mm3972/data/GRACE/GRACE_JPL_RL06/GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc"
    f          = addfile (filename,"r")

    tot_grace  = (GRACE_end_year-GRACE_start_year+1)*12 ; how many year will be read from GRACE
    GRACE(compensate_month+0:compensate_month+4,:,:)           = dble2flt(f->lwe_thickness(7:11,:,:))*10.
    GRACE(compensate_month+6:compensate_month+tot_grace-1,:,:) = dble2flt(f->lwe_thickness(12:tot_grace+5,:,:))*10.

    uncertainty(compensate_month+0:compensate_month+4,:,:)          = dble2flt(f->uncertainty(7:11,:,:))*10.
    uncertainty(compensate_month+6:compensate_month+tot_grace-1,:,:)= dble2flt(f->uncertainty(12:tot_grace+5,:,:))*10.

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

    res@tiMainString = "GRACE(compensate_month,:,:)"
    plot = gsn_csm_contour(wks,GRACE(compensate_month,:,:),res)

    do year = GRACE_start_year+1,GRACE_end_year
        ; print("add year :" + year)
        May_ave = May_ave + GRACE(compensate_month+(year-GRACE_start_year)*12+4,:,:)
        Jun_ave = Jun_ave + GRACE(compensate_month+(year-GRACE_start_year)*12+5,:,:)
        Jul_ave = Jul_ave + GRACE(compensate_month+(year-GRACE_start_year)*12+6,:,:)
    end do

    May_ave = May_ave/(GRACE_end_year-GRACE_start_year)
    Jun_ave = Jun_ave/(GRACE_end_year-GRACE_start_year)
    Jul_ave = Jul_ave/(GRACE_end_year-GRACE_start_year)
    ; print(May_ave)

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

end if
end if
end if

; ++++++++ check GRACE input ++++++++
res@tiMainString = "GRACE(compensate_month+24,:,:)"
plot = gsn_csm_contour(wks,GRACE(compensate_month+24,:,:),res)
res@tiMainString = "uncertainty(compensate_month+24,:,:)"
plot = gsn_csm_contour(wks,uncertainty(compensate_month+24,:,:),res)
; +++++++++++++++++++++++++++++++++++

; ====== post process GRACE =======
filename = "/g/data/w35/mm3972/model/cable/src/CABLE-AUX/offline/gridinfo_mmy_MD_elev_orig_std_avg-sand_landmask.nc"
f        = addfile (filename,"r")
landmask = where(f->landsea .eq. 0, 1, 0)

GRACE_lat = f->latitude
GRACE_lon = f->longitude

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

printVarSummary(GRACE_lat)
printVarSummary(GRACE_lon)
; regrid

GRACE_regrid       = linint2_Wrap(GRACE&lon,GRACE&lat,GRACE,True,CABLE_lon({C:D}),CABLE_lat({A:B}),0)
uncertainty_regrid = linint2_Wrap(uncertainty&lon,uncertainty&lat,uncertainty,True,CABLE_lon({C:D}),CABLE_lat({A:B}),0)
printVarSummary(GRACE_regrid)
printVarSummary(uncertainty_regrid)
; ++++++++ check GRACE input ++++++++
res@tiMainString = "dim_avg_n_Wrap(GRACE,0)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(GRACE,0),res)
res@tiMainString = "dim_avg_n_Wrap(GRACE_regrid,0)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(GRACE_regrid,0),res)
res@tiMainString = "dim_avg_n_Wrap(uncertainty,0)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(uncertainty,0),res)
res@tiMainString = "dim_avg_n_Wrap(uncertainty_regrid,0)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(uncertainty_regrid,0),res)
; +++++++++++++++++++++++++++++++++++

delete(GRACE)
delete(uncertainty)
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
          filename1 = path+case_name(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
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

; Rain_LIS_regrid   = regrid_data(Rain_LIS,interp_method, lat2d, lon2d, \
;                         lat, lon, src_mask, dst_mask )

; ++++++++ check GRACE input ++++++++
res@tiMainString = "dim_avg_n_Wrap(GRACE_regrid,0)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(GRACE_regrid,0),res)
res@tiMainString = "dim_avg_n_Wrap(TWS(0,:,{A:B},{C:D}),0),0)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(TWS(0,:,{A:B},{C:D}),0),res)
res@tiMainString = "dim_avg_n_Wrap(uncertainty_regrid,0),0)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(uncertainty_regrid,0),res)
; +++++++++++++++++++++++++++++++++++

; ======================= Calculate CABLE 2004-2009 avg ========================
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

; ; Calculate 2004-2009 TWS average
; if year_e .gt. 2009 then
;   adjust  = (year_e - GRACE_end_year)*12
; else
;   adjust  = 0
;   GRACE_end_year = year_e
; end if
tot_grace_mon = (2009-2004+1)*12
adjust        = (2004-year_s)*12
TWS_ave        =  dim_avg_n_Wrap(TWS(:,adjust:adjust+tot_grace_mon-1,:,:),1)
GWMoist_ave    =  dim_avg_n_Wrap(GWMoist(:,adjust:adjust+tot_grace_mon-1,:,:),1)
CanopInt_ave   =  dim_avg_n_Wrap(CanopInt(:,adjust:adjust+tot_grace_mon-1,:,:),1)
SWE_ave        =  dim_avg_n_Wrap(SWE(:,adjust:adjust+tot_grace_mon-1,:,:),1)
SoilMoist_ave  =  dim_avg_n_Wrap(SoilMoist(:,adjust:adjust+tot_grace_mon-1,:,:),1)

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
GRACE_global_ave       = new((/3,total/),float, -99999.)
uncertainty_global_ave = new(total,float, -99999.)
TWS_global_ave         = new((/case_sum,total/),float,default_value)
GWMoist_global_ave     = new((/case_sum,total/),float,default_value)
CanopInt_global_ave    = new((/case_sum,total/),float,default_value)
SWE_global_ave         = new((/case_sum,total/),float,default_value)
SoilMoist_global_ave   = new((/case_sum,total/),float,default_value)
Rain_global_ave        = new((/case_sum,total/),float,default_value)

GRACE_global_ave       = 0.0
uncertainty_global_ave = 0.0
TWS_global_ave         = 0.0
GWMoist_global_ave     = 0.0
CanopInt_global_ave    = 0.0
SWE_global_ave         = 0.0
SoilMoist_global_ave   = 0.0
Rain_global_ave        = 0.0

; print("average GRACE")
GRACE_tmp = where(PFT_mask_multi .eq. 1, GRACE_regrid,  -99999.)
GRACE_tmp@_FillValue = -99999.
uncertainty_tmp = where(PFT_mask_multi .eq. 1, uncertainty_regrid,  -99999.)
uncertainty_tmp@_FillValue = -99999.

GRACE_global_ave(0,:)     = dim_avg_n_Wrap(GRACE_tmp,(/1,2/))
uncertainty_global_ave    = dim_avg_n_Wrap(uncertainty_tmp,(/1,2/))

; ++++++++ check GRACE input ++++++++
res@tiMainString = "dim_avg_n_Wrap(GRACE_tmp,0),0)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(GRACE_tmp,0),res)
; +++++++++++++++++++++++++++++++++++
delete(GRACE_regrid)
delete(uncertainty_regrid)
delete(GRACE_tmp)
delete(uncertainty_tmp)

; print("average CABLE")
do case_num = 0,case_sum-1
  if cable_version .eq. "CABLE-HESS" then
      TWS_tmp = where(PFT_mask_multi .eq. 1, TWS(case_num,:,{A:B},{C:D}),  default_value)
      TWS_tmp@_FillValue = default_value
      GWMoist_tmp = where(PFT_mask_multi .eq. 1, GWMoist(case_num,:,{A:B},{C:D}),  default_value)
      GWMoist_tmp@_FillValue = default_value
      CanopInt_tmp = where(PFT_mask_multi .eq. 1, CanopInt(case_num,:,{A:B},{C:D}),  default_value)
      CanopInt_tmp@_FillValue = default_value
      SWE_tmp = where(PFT_mask_multi .eq. 1, SWE(case_num,:,{A:B},{C:D}),  default_value)
      SWE_tmp@_FillValue = default_value
      SoilMoist_tmp = where(PFT_mask_multi .eq. 1, SoilMoist(case_num,:,{A:B},{C:D}),  default_value)
      SoilMoist_tmp@_FillValue = default_value
      Rain_tmp = where(PFT_mask_multi .eq. 1, Rain(case_num,:,{A:B},{C:D}),  default_value)
      Rain_tmp@_FillValue = default_value

      TWS_global_ave(case_num,:)       = dim_avg_n_Wrap(TWS_tmp,(/1,2/))
      GWMoist_global_ave(case_num,:)   = dim_avg_n_Wrap(GWMoist_tmp,(/1,2/))
      CanopInt_global_ave(case_num,:)  = dim_avg_n_Wrap(CanopInt_tmp,(/1,2/))
      SWE_global_ave(case_num,:)       = dim_avg_n_Wrap(SWE_tmp,(/1,2/))
      SoilMoist_global_ave(case_num,:) = dim_avg_n_Wrap(SoilMoist_tmp,(/1,2/))
      Rain_global_ave(case_num,:)      = dim_avg_n_Wrap(Rain_tmp,(/1,2/))

      ; ++++++++ check GRACE input ++++++++
      res@tiMainString = "dim_avg_n_Wrap(TWS_tmp,0)"
      plot = gsn_csm_contour(wks,dim_avg_n_Wrap(TWS_tmp,0),res)
      ; +++++++++++++++++++++++++++++++++++

      delete(TWS_tmp)
      delete(GWMoist_tmp)
      delete(CanopInt_tmp)
      delete(SWE_tmp)
      delete(SoilMoist_tmp)
      delete(Rain_tmp)
  else if cable_version .eq. "LIS-CABLE" then
      TWS_global_ave(case_num,:)       = dim_avg_n_Wrap(TWS(case_num,time,:,:),(/1,2/))
      GWMoist_global_ave(case_num,:)   = dim_avg_n_Wrap(GWMoist(case_num,time,:,:),(/1,2/))
      CanopInt_global_ave(case_num,:)  = dim_avg_n_Wrap(CanopInt(case_num,time,:,:),(/1,2/))
      SWE_global_ave(case_num,:)       = dim_avg_n_Wrap(SWE(case_num,time,:,:),(/1,2/))
      SoilMoist_global_ave(case_num,:) = dim_avg_n_Wrap(SoilMoist(case_num,time,:,:),(/1,2/))
      Rain_global_ave(case_num,:)      = dim_avg_n_Wrap(Rain(case_num,time,:,:),(/1,2/))
  end if
  end if
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

delete(TWS)
delete(GWMoist)
delete(CanopInt)
delete(SWE)
delete(SoilMoist)
delete(Rain)

;GRACE & TWS is anormaly here
delete(wks)
delete(res)
delete(pic)

; __________________________ WRITE TXT OUTPUT _____________________________
option        = True
option@row    = False
; option@tspace = 0

do case_num = 0, case_sum-1
  option@fout = pic_message+"_GRACE_PFT-"+PFT+"_"+year_s+"-"+year_e+".txt"
  tmp      = new((/6,total/),float,default_value)
  tmp(0,:) = GRACE_global_ave(0,:)
  tmp(1,:) = TWS_global_ave(case_num,:)
  tmp(2,:) = GWMoist_global_ave(case_num,:)
  tmp(3,:) = SoilMoist_global_ave(case_num,:)
  tmp(4,:) = CanopInt_global_ave(case_num,:)
  tmp(5,:) = SWE_global_ave(case_num,:)
  print( tmp(0,0:11))
  printVarSummary( tmp)
  tmp = where (ismissing(tmp), -99999., tmp)
  tmp@_FillValue = -99999.
  print( tmp(0,0:11))
  format   = total+"f15.2"
  write_matrix (tmp, format, option)
  delete(tmp)
  delete(option@fout)
end do
end
