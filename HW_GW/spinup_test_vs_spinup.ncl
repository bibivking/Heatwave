;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

path = "/g/data/w35/mm3972/model/wrf/NUWRF/nuwrf_25km/LISWRF_configs/LIS_output"

var_name= (/\
            "Qle_tavg",\ ; "W m-2" "latent heat flux" ;
            "Qh_tavg",\  ; "W m-2" "sensible heat flux" ;
            "Qg_tavg",\  ; "W m-2" "soil heat flux" ;
            "Evap_tavg",\ ; "kg m-2 s-1" "total evapotranspiration" ;
            "TVeg_tavg",\ ; "kg m-2 s-1" "vegetation transpiration"
            "ESoil_tavg",\; "kg m-2 s-1" "bare soil evaporation"
            "Qs_tavg",\   ; "kg m-2 s-1" "surface runoff" ;
            "Qsb_tavg",\  ; "kg m-2 s-1" "subsurface runoff amount" ;
            "SnowDepth_inst",\; "m" "snow depth"
            "AvgSurfT_tavg",\; "K" "surface temperature" ;
            "SoilWet_inst",\ ; "-" "total soil wetness"
            "SoilTemp_inst",\;(time, 6, north_south, east_west) "K" "soil temperature" ;
            "SoilMoist_inst"\;(time, 6, north_south, east_west) "m^3 m-3" "soil moisture content"
           /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

;______________________________ Data _________________________________
year_s = 1979
year_e = 2008
year_sum = year_e-year_s+1
time_tot = 10958

Rain     = new((/time_tot/),float)
Qle      = new((/time_tot/),float)
Qh       = new((/time_tot/),float)
Qg       = new((/time_tot/),float)
Evap     = new((/time_tot/),float)
TVeg     = new((/time_tot/),float)
ESoil    = new((/time_tot/),float)
Qs       = new((/time_tot/),float)
Qsb      = new((/time_tot/),float)
SnowDepth= new((/time_tot/),float)
AvgSurfT = new((/time_tot/),float)
SoilWet  = new((/time_tot/),float)
SoilTemp = new((/time_tot,6/),float)
SoilMoist= new((/time_tot,6/),float)

i = 0

do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    do month = 1, 12
      print("month = "+month)
      if month .lt. 10 then
        filename = path+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename = path+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      print(filename)
      f = addfile (filename,"r")

      do day = 1, dom(month-1)
        print("day = "+day)

        Rain(i)     = avg(f->Rainf_f_inst(day-1,:,:))
        Qle(i)      = avg(f->Qle_tavg(day-1,:,:))
        Qh(i)       = avg(f->Qh_tavg(day-1,:,:))
        Qg(i)       = avg(f->Qg_tavg(day-1,:,:))
        Evap(i)     = avg(f->Evap_tavg(day-1,:,:))*3600*24
        TVeg(i)     = avg(f->TVeg_tavg(day-1,:,:))*3600*24
        ESoil(i)    = avg(f->ESoil_tavg(day-1,:,:))*3600*24
        Qs(i)       = avg(f->Qs_tavg(day-1,:,:))*3600*24
        Qsb(i)      = avg(f->Qsb_tavg(day-1,:,:))*3600*24
        SnowDepth(i)= avg(f->SnowDepth_inst(day-1,:,:))
        AvgSurfT(i) = avg(f->AvgSurfT_tavg(day-1,:,:))-273.15
        SoilWet(i)  = avg(f->SoilWet_inst(day-1,:,:))
        do soil = 0,5
          SoilTemp(i,soil) = avg(f->SoilTemp_inst(day-1,soil,:,:))-273.15
          SoilMoist(i,soil)= avg(f->SoilMoist_inst(day-1,soil,:,:))
        end do
        i = i + 1
      end do ; day
      delete(f)
    end do ; month
end do ; year

;_____________________________ PLOT ___________________________________
x = ispan(1,i,1)

pic = "Spinup_nuwrf_25km"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"BlueYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��
res                    = True
res@gsnMaximize        = True
res@vpHeightF          = 0.6
res@vpWidthF           = 0.6
res@vpXF               = 0.3
res@vpYF               = 0.85
res@pmLegendWidthF     = 0.14                     ; Resize legend width
res@pmLegendHeightF    = 0.11                     ; and height

res@xyLineThicknesses  = (/  1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0/)          ; make second line thicker
;res@pmLegendDisplayMode = "Always"
res@xyLineColors       = (/"blue","red","green","orange","red","red","red"/)          ; change line color
res@xyExplicitLegendLabels = (/"sim","trend","GW_1-60","NoAQ"/)
;res@trYMinF = 0

filename = path+"/LIS.CABLE.1999090100.d01.nc"
f = addfile(filename,"r")

Rain_rc    = regline(x, Rain)          ; slope
Rain_arry  = new ( (/2,time_tot/), typeof(Rain))
Rain_arry(0,:) = Rain
Rain_arry(1,:) = Rain_rc*x + Rain_rc@yintercept
res@tiYAxisString  = f->Rainf_f_inst@long_name+" (mm/d)"
res@tiXAxisString  = "days"
plot = gsn_csm_xy(wks,x,Rain_arry,res)

Qle_rc    = regline(x, Qle)          ; slope
Qle_arry  = new ( (/2,time_tot/), typeof(Qle))
Qle_arry(0,:) = Qle
Qle_arry(1,:) = Qle_rc*x + Qle_rc@yintercept
res@tiYAxisString  = f->Qle_tavg@long_name+" ("+f->Qle_tavg@units+")"
res@tiXAxisString  = "days"
plot = gsn_csm_xy(wks,x,Qle_arry,res)

Qh_rc    = regline(x, Qh)          ; slope
Qh_arry  = new ( (/2,time_tot/), typeof(Qh))
Qh_arry(0,:) = Qh
Qh_arry(1,:) = Qh_rc*x + Qh_rc@yintercept
res@tiYAxisString  = f->Qh_tavg@long_name+" ("+f->Qh_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qh_arry,res)

Qg_rc    = regline(x, Qg)          ; slope
Qg_arry  = new ( (/2,time_tot/), typeof(Qg))
Qg_arry(0,:) = Qg
Qg_arry(1,:) = Qg_rc*x + Qg_rc@yintercept
res@tiYAxisString  = f->Qg_tavg@long_name+" ("+f->Qg_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qg_arry,res)

Evap_rc    = regline(x, Evap)          ; slope
Evap_arry  = new ( (/2,time_tot/), typeof(Evap))
Evap_arry(0,:) = Evap
Evap_arry(1,:) = Evap_rc*x + Evap_rc@yintercept
res@tiYAxisString  = f->Evap_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Evap_arry,res)

TVeg_rc    = regline(x, TVeg)          ; slope
TVeg_arry  = new ( (/2,time_tot/), typeof(TVeg))
TVeg_arry(0,:) = TVeg
TVeg_arry(1,:) = TVeg_rc*x + TVeg_rc@yintercept
res@tiYAxisString  = f->TVeg_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,TVeg_arry,res)

ESoil_rc    = regline(x, ESoil)          ; slope
ESoil_arry  = new ( (/2,time_tot/), typeof(ESoil))
ESoil_arry(0,:) = ESoil
ESoil_arry(1,:) = ESoil_rc*x + ESoil_rc@yintercept
res@tiYAxisString  = f->ESoil_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,ESoil_arry,res)

Qs_rc    = regline(x, Qs)          ; slope
Qs_arry  = new ( (/2,time_tot/), typeof(Qs))
Qs_arry(0,:) = Qs
Qs_arry(1,:) = Qs_rc*x + Qs_rc@yintercept
res@tiYAxisString  = f->Qs_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qs_arry,res)

Qsb_rc    = regline(x, Qsb)          ; slope
Qsb_arry  = new ( (/2,time_tot/), typeof(Qsb))
Qsb_arry(0,:) = Qsb
Qsb_arry(1,:) = Qsb_rc*x + Qsb_rc@yintercept
res@tiYAxisString  = f->Qsb_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qsb_arry,res)

SnowDepth_rc    = regline(x, SnowDepth)          ; slope
SnowDepth_arry  = new ( (/2,time_tot/), typeof(SnowDepth))
SnowDepth_arry(0,:) = SnowDepth
SnowDepth_arry(1,:) = SnowDepth_rc*x + SnowDepth_rc@yintercept
res@tiYAxisString  = f->SnowDepth_inst@long_name+" ("+f->SnowDepth_inst@units+")"
plot = gsn_csm_xy(wks,x,SnowDepth_arry,res)

do soil = 0,5
  res@trYMaxF = 50.
  res@trYMinF = 0.
  SoilTemp_rc    = regline(x, SoilTemp(:,soil))          ; slope
  SoilTemp_arry  = new ( (/2,time_tot/), typeof(SoilTemp))
  SoilTemp_arry(0,:) = SoilTemp(:,soil)
  SoilTemp_arry(1,:) = SoilTemp_rc*x + SoilTemp_rc@yintercept
  res@tiYAxisString  = "Layer "+(soil+1)+" Soil Temp (C)"
  plot = gsn_csm_xy(wks,x,SoilTemp_arry,res)
  delete(SoilTemp_rc)
  delete(SoilTemp_arry)
  delete(res@trYMaxF)
  delete(res@trYMinF)
end do
do soil = 0,5
  res@trYMaxF = 0.4
  res@trYMinF = 0.
  SoilMoist_rc    = regline(x, SoilMoist(:,soil))          ; slope
  SoilMoist_arry  = new ( (/2,time_tot/), typeof(SoilMoist))
  SoilMoist_arry(0,:) = SoilMoist(:,soil)
  SoilMoist_arry(1,:) = SoilMoist_rc*x + SoilMoist_rc@yintercept
  res@tiYAxisString  = "Layer "+(soil+1)+" Soil Moisture (m3/m3)"
  plot = gsn_csm_xy(wks,x,SoilMoist_arry,res)
  delete(SoilMoist_rc)
  delete(SoilMoist_arry)
  delete(res@trYMaxF)
  delete(res@trYMinF)
end do

res@trYMaxF = 50.
res@trYMinF = 0.
AvgSurfT_rc    = regline(x, AvgSurfT)          ; slope
AvgSurfT_arry  = new ( (/2,time_tot/), typeof(AvgSurfT))
AvgSurfT_arry(0,:) = AvgSurfT
AvgSurfT_arry(1,:) = AvgSurfT_rc*x + AvgSurfT_rc@yintercept
res@tiYAxisString  = f->AvgSurfT_tavg@long_name+" (C)"
plot = gsn_csm_xy(wks,x,AvgSurfT_arry,res)
delete(res@trYMaxF)
delete(res@trYMinF)

SoilWet_rc    = regline(x, SoilWet)          ; slope
SoilWet_arry  = new ( (/2,time_tot/), typeof(SoilWet))
SoilWet_arry(0,:) = SoilWet
SoilWet_arry(1,:) = SoilWet_rc*x + SoilWet_rc@yintercept
res@tiYAxisString = f->SoilWet_inst@long_name+" ("+f->SoilWet_inst@units+")"
plot = gsn_csm_xy(wks,x,SoilWet_arry,res)

end
