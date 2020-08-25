;*****************************************
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


;______________________________ Data _________________________________
path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/Princeton_ctl/LIS_output"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/ctl_para/outputs"

year_s = 2003
year_e = 2008
year_sum = year_e-year_s+1
day_sum  = 6940

Soil_thickness  = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Soil_thickness1 = Soil_thickness
;(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE

Var     = new((/2,19/),float)
Var     = 0.0

Var_name = (/\
             "Rain","Evap","TVeg","ESoil","ECanop","Runoff",\
             "Rnet","Qle","Qh","Qg","Tair","SoilMoist","SoilTemp",\
             "SoilMoist1","SoilMoist2","SoilMoist3","SoilMoist4", \
             "SoilMoist5","SoilMoist6"\
             /)

; ================================ Set Regrid ==================================
regrid = False
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

; ================================= READ DATA ==================================

if regrid then

  ;  READ CABLE
  filename = path1+"/LIS.CABLE."+year_s+"010100.d01.nc"
  f        = addfile (filename,"r")

  Rain_LIS     = dim_sum_n_Wrap(f->Rainf_tavg,0) ;
  Tair_LIS     = dim_avg_n_Wrap(f->Tair_f_tavg,0)
  Qair_LIS     = dim_avg_n_Wrap(f->Qair_f_tavg,0)
  Psurf_LIS    = dim_avg_n_Wrap(f->Psurf_f_tavg,0)
  SWdown_LIS   = dim_avg_n_Wrap(f->SWdown_f_tavg,0)
  LWdown_LIS   = dim_avg_n_Wrap(f->LWdown_f_tavg,0)
  Wind_LIS     = dim_avg_n_Wrap(f->Wind_f_tavg,0)
  Evap_LIS     = dim_sum_n_Wrap(f->Evap_tavg,0)
  TVeg_LIS     = dim_sum_n_Wrap(f->TVeg_tavg,0)
  ESoil_LIS    = dim_sum_n_Wrap(f->ESoil_tavg,0)
  ECanop_LIS   = dim_sum_n_Wrap(f->ECanop_tavg,0)
  Runoff_LIS   = dim_sum_n_Wrap(f->Qs_tavg,0)
  Rnet_LIS     = dim_avg_n_Wrap(f->Swnet_tavg,0)
  Qle_LIS      = dim_avg_n_Wrap(f->Qle_tavg,0)
  Qh_LIS       = dim_avg_n_Wrap(f->Qh_tavg,0)
  Qg_LIS       = dim_avg_n_Wrap(f->Qg_tavg,0)
  SMsurf_LIS   = dim_avg_n_Wrap(f->SoilMoist_tavg(:,0,:,:),0)
  SoilTemp_LIS = dim_avg_n_Wrap(f->SoilTemp_tavg(:,0,:,:),0)
  SoilMoist_LIS= dim_avg_n_Wrap(f->SoilMoist_tavg(:,0,:,:),0)
  LAI_LIS      = dim_avg_n_Wrap(f->LAI_inst,0)

  Rain_LIS     = 0.
  Tair_LIS     = 0.
  Qair_LIS     = 0.
  Psurf_LIS    = 0.
  SWdown_LIS   = 0.
  LWdown_LIS   = 0.
  Wind_LIS     = 0.
  Evap_LIS     = 0.
  TVeg_LIS     = 0.
  ESoil_LIS    = 0.
  ECanop_LIS   = 0.
  Runoff_LIS   = 0.
  Rnet_LIS     = 0.
  Qle_LIS      = 0.
  Qh_LIS       = 0.
  Qg_LIS       = 0.
  SoilTemp_LIS = 0.
  SoilMoist_LIS= 0.
  SMsurf_LIS   = 0.
  LAI_LIS      = 0.
  delete(f)

  filename1 = path2+"/cable_out_"+year_s+".nc"
  f1        = addfile (filename1,"r")

  Rain_HESS     = dim_sum_n_Wrap(f1->Rainf(:,{A:B},{C:D}),0)
  Tair_HESS     = dim_avg_n_Wrap(f1->Tair(:,{A:B},{C:D}),0)
  Qair_HESS     = dim_avg_n_Wrap(f1->Qair(:,{A:B},{C:D}),0)
  Psurf_HESS    = dim_avg_n_Wrap(f1->PSurf(:,{A:B},{C:D}),0)
  SWdown_HESS   = dim_avg_n_Wrap(f1->SWdown(:,{A:B},{C:D}),0)
  LWdown_HESS   = dim_avg_n_Wrap(f1->LWdown(:,{A:B},{C:D}),0)
  Wind_HESS     = dim_avg_n_Wrap(f1->Wind(:,{A:B},{C:D}),0)

  Evap_HESS     = dim_sum_n_Wrap(f1->Evap(:,{A:B},{C:D}),0)
  TVeg_HESS     = dim_sum_n_Wrap(f1->TVeg(:,{A:B},{C:D}),0)
  ESoil_HESS    = dim_sum_n_Wrap(f1->ESoil(:,{A:B},{C:D}),0)
  ECanop_HESS   = dim_sum_n_Wrap(f1->ECanop(:,{A:B},{C:D}),0)
  Runoff_HESS   = dim_sum_n_Wrap(f1->Qs(:,{A:B},{C:D}),0)
  Rnet_HESS     = dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D}),0)
  Qle_HESS      = dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D}),0)
  Qh_HESS       = dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D}),0)
  Qg_HESS       = dim_avg_n_Wrap(f1->Qg(:,{A:B},{C:D}),0)
  SMsurf_HESS   = dim_avg_n_Wrap(f1->SoilMoist(:,0,{A:B},{C:D}),0)
  SoilTemp_HESS = dim_avg_n_Wrap(f1->SoilTemp(:,0,{A:B},{C:D}),0)
  SoilMoist_HESS= dim_avg_n_Wrap(f1->SoilMoist(:,0,{A:B},{C:D}),0)
  LAI_HESS      = dim_avg_n_Wrap(f1->LAI(:,{A:B},{C:D}),0)

  Rain_HESS     = 0.
  Tair_HESS     = 0.
  Qair_HESS     = 0.
  Psurf_HESS    = 0.
  SWdown_HESS   = 0.
  LWdown_HESS   = 0.
  Wind_HESS     = 0.
  Evap_HESS     = 0.
  TVeg_HESS     = 0.
  ESoil_HESS    = 0.
  ECanop_HESS   = 0.
  Runoff_HESS   = 0.
  Rnet_HESS     = 0.
  Qle_HESS      = 0.
  Qh_HESS       = 0.
  Qg_HESS       = 0.
  SoilTemp_HESS = 0.
  SoilMoist_HESS= 0.
  SMsurf_HESS   = 0.
  LAI_HESS      = 0.
  delete(f1)

  do year = year_s,year_e
      print("year = " +year)

      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if

      filename1 = path2+"/cable_out_"+year+".nc"
      f1 = addfile (filename1,"r")

      ; _________________________ read LIS-CABLE ____________________________
      do month = 1, 12
        if month .lt. 10 then
          filename = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
        else
          filename = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
        end if
        f = addfile (filename,"r")

        Rain_LIS     = Rain_LIS   + dim_sum_n_Wrap(f->Rainf_tavg,0)*3600*24
        Evap_LIS     = Evap_LIS   + dim_sum_n_Wrap(f->Evap_tavg,0)*3600*24
        TVeg_LIS     = TVeg_LIS   + dim_sum_n_Wrap(f->TVeg_tavg,0)*3600*24
        ESoil_LIS    = ESoil_LIS  + dim_sum_n_Wrap(f->ESoil_tavg,0)*3600*24
        ECanop_LIS   = ECanop_LIS + dim_sum_n_Wrap(f->ECanop_tavg,0)*3600*24
        Runoff_LIS   = Runoff_LIS + dim_sum_n_Wrap(f->Qs_tavg,0)*3600*24 \
                                  + dim_sum_n_Wrap(f->Qsb_tavg,0)*3600*24
        Tair_LIS     = Tair_LIS   + dim_avg_n_Wrap(f->Tair_f_tavg,0)-273.15
        Qair_LIS     = Qair_LIS   + dim_avg_n_Wrap(f->Qair_f_tavg,0)
        Psurf_LIS    = Psurf_LIS  + dim_avg_n_Wrap(f->Psurf_f_tavg,0)/100.
        SWdown_LIS   = SWdown_LIS + dim_avg_n_Wrap(f->SWdown_f_tavg,0)
        LWdown_LIS   = LWdown_LIS + dim_avg_n_Wrap(f->LWdown_f_tavg,0)
        Wind_LIS     = Wind_LIS   + dim_avg_n_Wrap(f->Wind_f_tavg,0)
        Rnet_LIS     = Rnet_LIS   + dim_avg_n_Wrap(f->Swnet_tavg,0)\
                                  + dim_avg_n_Wrap(f->Lwnet_tavg,0)
        Qle_LIS      = Qle_LIS    + dim_avg_n_Wrap(f->Qle_tavg,0)
        Qh_LIS       = Qh_LIS     + dim_avg_n_Wrap(f->Qh_tavg,0)
        Qg_LIS       = Qg_LIS     + dim_avg_n_Wrap(f->Qg_tavg,0)
        LAI_LIS      = LAI_LIS    + dim_avg_n_Wrap(f->LAI_inst,0)
        SMsurf_LIS   = SMsurf_LIS + ( dim_avg_n_Wrap(f->SoilMoist_tavg(:,0,:,:),0)*0.005 \
                                    + dim_avg_n_Wrap(f->SoilMoist_tavg(:,1,:,:),0)*0.075 \
                                    + dim_avg_n_Wrap(f->SoilMoist_tavg(:,2,:,:),0)*(0.1 - 0.005 - 0.075)) \
                                    / 0.1
        do soil = 0,5
          SoilTemp_LIS = SoilTemp_LIS + (dim_avg_n_Wrap(f->SoilTemp_tavg(:,soil,:,:),0)-273.15)\
                                      * Soil_thickness1(soil)/sum(Soil_thickness1)
          SoilMoist_LIS= SoilMoist_LIS +dim_avg_n_Wrap(f->SoilMoist_tavg(:,soil,:,:),0)\
                                         *Soil_thickness1(soil)/sum(Soil_thickness1)
        end do ; soil
        delete(f)

        ; _________________________ read CABLE (HESS) _________________________

        Rain_HESS     = Rain_HESS   + f1->Rainf(month-1,{A:B},{C:D})*dom(month-1)*3600*24
        Evap_HESS     = Evap_HESS   + f1->Evap(month-1,{A:B},{C:D})*dom(month-1)*3600*24
        TVeg_HESS     = TVeg_HESS   + f1->TVeg(month-1,{A:B},{C:D})*dom(month-1)*3600*24
        ESoil_HESS    = ESoil_HESS  + f1->ESoil(month-1,{A:B},{C:D})*dom(month-1)*3600*24
        ECanop_HESS   = ECanop_HESS + f1->ECanop(month-1,{A:B},{C:D})*dom(month-1)*3600*24
        Runoff_HESS   = Runoff_HESS + f1->Qs(month-1,{A:B},{C:D})*dom(month-1)*3600*24 \
                                    + f1->Qsb(month-1,{A:B},{C:D})*dom(month-1)*3600*24

      end do ; month

      Tair_HESS     = Tair_HESS   + dim_avg_n_Wrap(f1->Tair(:,{A:B},{C:D}),0)-273.15
      Qair_HESS     = Qair_HESS   + dim_avg_n_Wrap(f1->Qair(:,{A:B},{C:D}),0)
      Psurf_HESS    = Psurf_HESS  + dim_avg_n_Wrap(f1->PSurf(:,{A:B},{C:D}),0)
      SWdown_HESS   = SWdown_HESS + dim_avg_n_Wrap(f1->SWdown(:,{A:B},{C:D}),0)
      LWdown_HESS   = LWdown_HESS + dim_avg_n_Wrap(f1->LWdown(:,{A:B},{C:D}),0)
      Wind_HESS     = Wind_HESS   + dim_avg_n_Wrap(f1->Wind(:,{A:B},{C:D}),0)
      Rnet_HESS     = Rnet_HESS   + dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D}),0)
      Qle_HESS      = Qle_HESS    + dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D}),0)
      Qh_HESS       = Qh_HESS     + dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D}),0)
      Qg_HESS       = Qg_HESS     + dim_avg_n_Wrap(f1->Qg(:,{A:B},{C:D}),0)
      LAI_HESS      = LAI_HESS    + dim_avg_n_Wrap(f1->LAI(:,{A:B},{C:D}),0)
      SMsurf_HESS   = SMsurf_HESS + ( dim_avg_n_Wrap(f1->SoilMoist(:,0,{A:B},{C:D}),0)*0.022 \
                                  + dim_avg_n_Wrap(f1->SoilMoist(:,1,{A:B},{C:D}),0)*0.058 \
                                  + dim_avg_n_Wrap(f1->SoilMoist(:,2,{A:B},{C:D}),0)*(0.1 - 0.022 - 0.058)) \
                                  / 0.1
      do soil = 0,5
        SoilTemp_HESS = SoilTemp_HESS + ( dim_avg_n_Wrap(f1->SoilTemp(:,soil,{A:B},{C:D}),0)-273.15)\
                                   *Soil_thickness(soil)/sum(Soil_thickness)
        SoilMoist_HESS= SoilMoist_HESS + ( dim_avg_n_Wrap(f1->SoilMoist(:,soil,{A:B},{C:D}),0))\
                                   *Soil_thickness(soil)/sum(Soil_thickness)
      end do ; soil
      delete(f1)
  end do ; year
  Rain_LIS     = Rain_LIS/year_sum
  Evap_LIS     = Evap_LIS/year_sum
  TVeg_LIS     = TVeg_LIS/year_sum
  ESoil_LIS    = ESoil_LIS/year_sum
  ECanop_LIS   = ECanop_LIS/year_sum
  Runoff_LIS   = Runoff_LIS/year_sum
  Tair_LIS     = Tair_LIS/year_sum/12
  Qair_LIS     = Qair_LIS/year_sum/12
  Psurf_LIS    = Psurf_LIS/year_sum/12
  SWdown_LIS   = SWdown_LIS/year_sum/12
  LWdown_LIS   = LWdown_LIS/year_sum/12
  Wind_LIS     = Wind_LIS/year_sum/12
  Rnet_LIS     = Rnet_LIS/year_sum/12
  Qle_LIS      = Qle_LIS/year_sum/12
  Qh_LIS       = Qh_LIS/year_sum/12
  Qg_LIS       = Qg_LIS/year_sum/12
  SMsurf_LIS   = SMsurf_LIS/year_sum/12
  LAI_LIS      = LAI_LIS/year_sum/12
  SoilTemp_LIS = SoilTemp_LIS/year_sum/12
  SoilMoist_LIS= SoilMoist_LIS/year_sum/12

  Rain_HESS     = Rain_HESS/year_sum
  Evap_HESS     = Evap_HESS/year_sum
  TVeg_HESS     = TVeg_HESS/year_sum
  ESoil_HESS    = ESoil_HESS/year_sum
  ECanop_HESS   = ECanop_HESS/year_sum
  Runoff_HESS   = Runoff_HESS/year_sum
  Tair_HESS     = Tair_HESS/year_sum
  Qair_HESS     = Qair_HESS/year_sum
  Psurf_HESS    = Psurf_HESS/year_sum
  SWdown_HESS   = SWdown_HESS/year_sum
  LWdown_HESS   = LWdown_HESS/year_sum
  Wind_HESS     = Wind_HESS/year_sum
  Rnet_HESS     = Rnet_HESS/year_sum
  Qle_HESS      = Qle_HESS/year_sum
  Qh_HESS       = Qh_HESS/year_sum
  Qg_HESS       = Qg_HESS/year_sum
  LAI_HESS      = LAI_HESS/year_sum
  SMsurf_HESS   = SMsurf_HESS/year_sum
  SoilTemp_HESS = SoilTemp_HESS/year_sum
  SoilMoist_HESS= SoilMoist_HESS/year_sum

  ; ================================== Regrid ====================================
  Rain_LIS_regrid   = regrid_data(Rain_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Evap_LIS_regrid   = regrid_data(Evap_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  TVeg_LIS_regrid   = regrid_data(TVeg_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  ESoil_LIS_regrid  = regrid_data(ESoil_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  ECanop_LIS_regrid = regrid_data(ECanop_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Runoff_LIS_regrid = regrid_data(Runoff_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Tair_LIS_regrid   = regrid_data(Tair_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Qair_LIS_regrid   = regrid_data(Qair_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Psurf_LIS_regrid  = regrid_data(Psurf_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  SWdown_LIS_regrid = regrid_data(SWdown_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  LWdown_LIS_regrid = regrid_data(LWdown_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Wind_LIS_regrid   = regrid_data(Wind_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Rnet_LIS_regrid   = regrid_data(Rnet_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Qle_LIS_regrid    = regrid_data(Qle_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Qh_LIS_regrid     = regrid_data(Qh_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  Qg_LIS_regrid     = regrid_data(Qg_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  SMsurf_LIS_regrid = regrid_data(SMsurf_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  LAI_LIS_regrid    = regrid_data(LAI_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  SoilTemp_LIS_regrid  = regrid_data(SoilTemp_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )
  SoilMoist_LIS_regrid = regrid_data(SoilMoist_LIS,interp_method, lat2d, lon2d, \
                          lat, lon, src_mask, dst_mask )

  ; 0: Rain
  Var(0,0)     = avg(Rain_LIS_regrid)
  ; 1: Evap
  Var(0,1)     = avg(Evap_LIS_regrid)
  ; 2: TVeg
  Var(0,2)     = avg(TVeg_LIS_regrid)
  ; 3: ESoil
  Var(0,3)     = avg(ESoil_LIS_regrid)
  ; 4: ECanop
  Var(0,4)     = avg(ECanop_LIS_regrid)
  ; 5: runoff
  Var(0,5)     = avg(Runoff_LIS_regrid)
  ; 6: Rnet
  Var(0,6)     = avg(Rnet_LIS_regrid)
  ; 7: Qle
  Var(0,7)     = avg(Qle_LIS_regrid)
  ; 8: Qh
  Var(0,8)     = avg(Qh_LIS_regrid)
  ; 9: Qg
  Var(0,9)     = avg(Qg_LIS_regrid)
  ; 10: Tair
  Var(0,10)    = avg(Tair_LIS_regrid)
  ; 11: SoilMoist
  Var(0,11)     = avg(SoilMoist_LIS_regrid)
  ; 12: SoilTemp
  Var(0,12)     = avg(SoilTemp_LIS_regrid)

  do soil = 0,5
    ; 13-18: SoilMoist1-6
    Var(0,13+soil)= avg(SoilMoist_LIS_regrid)
  end do ; soil

  ; 0: Rain
  Var(1,0)     = avg(Rain_HESS)
  ; 1: Evap
  Var(1,1)     = avg(Evap_HESS)
  ; 2: TVeg
  Var(1,2)     = avg(TVeg_HESS)
  ; 3: ESoil
  Var(1,3)     = avg(ESoil_HESS)
  ; 4: ECanop
  Var(1,4)     = avg(ECanop_HESS)
  ; 5: runoff
  Var(1,5)     = avg(Runoff_HESS)
  ; 6: Rnet
  Var(1,6)     = avg(Rnet_HESS)
  ; 7: Qle
  Var(1,7)     = avg(Qle_HESS)
  ; 8: Qh
  Var(1,8)     = avg(Qh_HESS)
  ; 9: Qg
  Var(1,9)     = avg(Qg_HESS)
  ; 10: Tair
  Var(1,10)    = avg(Tair_HESS)
  ; 11: SoilMoist
  Var(1,11)    = avg(SoilMoist_HESS)
  ; 12: SoilTemp
  Var(1,12)    = avg(SoilTemp_HESS)

  do soil = 0,5
    ; 13-18: SoilMoist1-6
    Var(1,13+soil)= avg(SoilMoist_HESS)
  end do ; soil
else
  do year = year_s,year_e

      print("year = " +year)
      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if

      ; read CABLE (HESS)
      filename1 = path2+"/cable_out_"+year+".nc"
      f1 = addfile (filename1,"r")

      ; read LIS-CABLE
      do month = 1, 12
        ;print("month = "+month)
        if month .lt. 10 then
          filename = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
        else
          filename = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
        end if

        ; print(filename)
        f = addfile (filename,"r")

        ; 0: Rain
        Var(0,0)     = Var(0,0) + avg(dim_sum_n_Wrap(f->Rainf_tavg,0))*3600*24
        ; 1: Evap
        Var(0,1)     = Var(0,1) + avg(dim_sum_n_Wrap(f->Evap_tavg,0))*3600*24
        ; 2: TVeg
        Var(0,2)     = Var(0,2) + avg(dim_sum_n_Wrap(f->TVeg_tavg,0))*3600*24
        ; 3: ESoil
        Var(0,3)     = Var(0,3) + avg(dim_sum_n_Wrap(f->ESoil_tavg,0))*3600*24
        ; 4: ECanop
        Var(0,4)     = Var(0,4) + avg(dim_sum_n_Wrap(f->ECanop_tavg,0))*3600*24
        ; 5: Qs
        Var(0,5)     = Var(0,5) + avg(dim_sum_n_Wrap(f->Qs_tavg,0))*3600*24\
                                + avg(dim_sum_n_Wrap(f->Qsb_tavg,0))*3600*24
        ; 6: Rnet
        Var(0,6)     = Var(0,6) + avg(dim_avg_n_Wrap(f->Swnet_tavg,0))\
                                + avg(dim_avg_n_Wrap(f->Lwnet_tavg,0))
        ; 7: Qle
        Var(0,7)     = Var(0,7) + avg(dim_avg_n_Wrap(f->Qle_tavg,0))
        ; 8: Qh
        Var(0,8)     = Var(0,8) + avg(dim_avg_n_Wrap(f->Qh_tavg,0))
        ; 9: Qg
        Var(0,9)     = Var(0,9) + avg(dim_avg_n_Wrap(f->Qg_tavg,0))
        ; 10: Tair
        Var(0,10)     = Var(0,10) + avg(dim_avg_n_Wrap(f->Tair_f_tavg,0))-273.15

        do soil = 0,5
          ; 11: SoilMoist
          Var(0,11)     = Var(0,11) + avg(dim_avg_n_Wrap(f->SoilMoist_tavg(:,soil,:,:),0))\
                                     *Soil_thickness1(soil)/sum(Soil_thickness1)
          ; 12: SoilTemp
          Var(0,12)     = Var(0,12) + (avg(dim_avg_n_Wrap(f->SoilTemp_tavg(:,soil,:,:),0))-273.15)\
                                     *Soil_thickness1(soil)/sum(Soil_thickness1)
          ; 13-18: SoilMoist1-6
          Var(0,13+soil)= Var(0,13+soil) + avg(dim_avg_n_Wrap(f->SoilMoist_tavg(:,soil,:,:),0))
        end do ; soil
        delete(f)


        ; read CABLE (HESS)

        ; 0: Rain
        Var(1,0)     = Var(1,0) + avg(f1->Rainf(month-1,{A:B},{C:D})*dom(month-1))*3600*24
        ; 1: Evap
        Var(1,1)     = Var(1,1) + avg(f1->Evap(month-1,{A:B},{C:D})*dom(month-1))*3600*24
        ; 2: TVeg
        Var(1,2)     = Var(1,2) + avg(f1->TVeg(month-1,{A:B},{C:D})*dom(month-1))*3600*24
        ; 3: ESoil
        Var(1,3)     = Var(1,3) + avg(f1->ESoil(month-1,{A:B},{C:D})*dom(month-1))*3600*24
        ; 4: ECanop
        Var(1,4)     = Var(1,4) + avg(f1->ECanop(month-1,{A:B},{C:D})*dom(month-1))*3600*24
        ; 5: Qs
        Var(1,5)     = Var(1,5) + avg(f1->Qs(month-1,{A:B},{C:D})*dom(month-1))*3600*24 \
                                + avg(f1->Qsb(month-1,{A:B},{C:D})*dom(month-1))*3600*24


      end do ; month

      ; 6: Rnet
      Var(1,6)     = Var(1,6) + avg(dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D}),0))
      ; 7: Qle
      Var(1,7)     = Var(1,7) + avg(dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D}),0))
      ; 8: Qh
      Var(1,8)     = Var(1,8) + avg(dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D}),0))
      ; 9: Qg
      Var(1,9)     = Var(1,9) + avg(dim_avg_n_Wrap(f1->Qg(:,{A:B},{C:D}),0))
      ; 10: Tair
      Var(1,10)     = Var(1,10) + avg(dim_avg_n_Wrap(f1->Tair(:,{A:B},{C:D}),0))-273.15

      do soil = 0,5
        ; 11: SoilMoist
        Var(1,11)     = Var(1,11) + avg(dim_avg_n_Wrap(f1->SoilMoist(:,soil,{A:B},{C:D}),0))\
                                   *Soil_thickness(soil)/sum(Soil_thickness)
        ; 12: SoilTemp
        Var(1,12)     = Var(1,12) + (avg(dim_avg_n_Wrap(f1->SoilTemp(:,soil,{A:B},{C:D}),0))-273.15)\
                                   *Soil_thickness(soil)/sum(Soil_thickness)
        ; 13-18: SoilMoist1-6
        Var(1,13+soil)= Var(1,13+soil) + avg(dim_avg_n_Wrap(f1->SoilMoist(:,soil,{A:B},{C:D}),0))
      end do ; soil
      delete(f1)
  end do ; year
  print(Var)
  Var(0,0:5) = Var(0,0:5)/year_sum
  Var(0,6:) = Var(0,6:)/year_sum/12
  Var(1,0:5) = Var(1,0:5)/year_sum
  Var(1,6:) = Var(1,6:)/year_sum
end if

do var_num = 0, dimsizes(Var_name)-1
  print(Var_name(var_num)+ " "+ Var(0,var_num)+" "+ Var(1,var_num))
end do

end
