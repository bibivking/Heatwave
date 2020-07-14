;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 500000000000 ;
end setvalues
;_____________________________________


begin

path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/nuwrf_25km/LISWRF_configs/LIS_output"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/outputs"


Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers


; region

A = -43.4377
B = -13.5853
C = 117.552
D = 153.598

;______________________________ Data _________________________________
year_s = 1979
year_e = 2008
year_sum = year_e-year_s+1
time_tot = 10958

Rain     = new((/2,time_tot/),float)
Qle      = new((/2,time_tot/),float)
Qh       = new((/2,time_tot/),float)
Qg       = new((/2,time_tot/),float)
Evap     = new((/2,time_tot/),float)
TVeg     = new((/2,time_tot/),float)
ESoil    = new((/2,time_tot/),float)
Qs       = new((/2,time_tot/),float)
Qsb      = new((/2,time_tot/),float)
SnowDepth= new((/2,time_tot/),float)
AvgSurfT = new((/2,time_tot/),float)
SoilTemp = new((/2,time_tot,6/),float)
SoilMoist= new((/2,time_tot,6/),float)

i = 0
j = 0
do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    ; read LIS-CABLE
    do month = 1, 12
      print("month = "+month)
      if month .lt. 10 then
        filename = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      ;print(filename)
      f = addfile (filename,"r")

      do day = 1, dom(month-1)
        ; print("day = "+day)

        Rain(0,i)     = avg(f->Rainf_f_inst(day-1,:,:))
        Qle(0,i)      = avg(f->Qle_tavg(day-1,:,:))
        Qh(0,i)       = avg(f->Qh_tavg(day-1,:,:))
        Qg(0,i)       = avg(f->Qg_tavg(day-1,:,:))
        Evap(0,i)     = avg(f->Evap_tavg(day-1,:,:))*3600*24
        TVeg(0,i)     = avg(f->TVeg_tavg(day-1,:,:))*3600*24
        ESoil(0,i)    = avg(f->ESoil_tavg(day-1,:,:))*3600*24
        Qs(0,i)       = avg(f->Qs_tavg(day-1,:,:))*3600*24
        Qsb(0,i)      = avg(f->Qsb_tavg(day-1,:,:))*3600*24
        SnowDepth(0,i)= avg(f->SnowDepth_inst(day-1,:,:))
        AvgSurfT(0,i) = avg(f->AvgSurfT_tavg(day-1,:,:))-273.15
        do soil = 0,5
          SoilTemp(0,i,soil) = avg(f->SoilTemp_inst(day-1,soil,:,:))-273.15
          SoilMoist(0,i,soil)= avg(f->SoilMoist_inst(day-1,soil,:,:))
        end do
        i = i + 1
      end do ; day
      delete(filename)
      delete(f)
    end do ; month

    ; read CABLE (HESS)
    filename1 = path2+"/cable_out_"+year+".nc"
    f1 = addfile (filename1,"r")
    do day = 1,dimsizes(f1->time)
      Rain(1,j)     = avg(f1->Rainf(day-1,{A:B},{C:D}))
      Qle(1,j)      = avg(f1->Qle(day-1,{A:B},{C:D}))
      Qh(1,j)       = avg(f1->Qh(day-1,{A:B},{C:D}))
      Qg(1,j)       = avg(f1->Qg(day-1,{A:B},{C:D}))
      Evap(1,j)     = avg(f1->Evap(day-1,{A:B},{C:D}))*3600*24
      TVeg(1,j)     = avg(f1->TVeg(day-1,{A:B},{C:D}))*3600*24
      ESoil(1,j)    = avg(f1->ESoil(day-1,{A:B},{C:D}))*3600*24
      Qs(1,j)       = avg(f1->Qs(day-1,{A:B},{C:D}))*3600*24
      Qsb(1,j)      = avg(f1->Qsb(day-1,{A:B},{C:D}))*3600*24
      SnowDepth(1,j)= avg(f1->SnowDepth(day-1,{A:B},{C:D}))
      AvgSurfT(1,j) = avg(f1->Tair(day-1,{A:B},{C:D}))-273.15
      do soil = 0,5
        SoilTemp(1,j,soil) = avg(f1->SoilTemp(day-1,soil,{A:B},{C:D}))-273.15
        SoilMoist(1,j,soil)= avg(f1->SoilMoist(day-1,soil,{A:B},{C:D}))
      end do
      j = j + 1
    end do ; day
    delete(f1)
end do ; year

;_____________________________ PLOT ___________________________________
x = ispan(1,i,1)

pic = "Spinup_LIS-CABLE_vs_CABLE-HESS"
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
res@xyExplicitLegendLabels = (/"LIS-CABLE","CABLE-HESS","GW_1-60","NoAQ"/)
;res@trYMinF = 0

filename = path1+"/LIS.CABLE.1999090100.d01.nc"
f = addfile(filename,"r")

; Rain_rc    = regline(x, Rain)          ; slope
; Rain_arry  = new ( (/2,time_tot/), typeof(Rain))
; Rain_arry(0,:) = Rain
; Rain_arry(1,:) = Rain_rc*x + Rain_rc@yintercept
res@tiYAxisString  = f->Rainf_f_inst@long_name+" (mm/d)"
res@tiXAxisString  = "days"
plot = gsn_csm_xy(wks,x,Rain,res)

; Qle_rc    = regline(x, Qle)          ; slope
; Qle_arry  = new ( (/2,time_tot/), typeof(Qle))
; Qle_arry(0,:) = Qle
; Qle_arry(1,:) = Qle_rc*x + Qle_rc@yintercept
res@tiYAxisString  = f->Qle_tavg@long_name+" ("+f->Qle_tavg@units+")"
res@tiXAxisString  = "days"
plot = gsn_csm_xy(wks,x,Qle,res)

; Qh_rc    = regline(x, Qh)          ; slope
; Qh_arry  = new ( (/2,time_tot/), typeof(Qh))
; Qh_arry(0,:) = Qh
; Qh_arry(1,:) = Qh_rc*x + Qh_rc@yintercept
res@tiYAxisString  = f->Qh_tavg@long_name+" ("+f->Qh_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qh,res)

; Qg_rc    = regline(x, Qg)          ; slope
; Qg_arry  = new ( (/2,time_tot/), typeof(Qg))
; Qg_arry(0,:) = Qg
; Qg_arry(1,:) = Qg_rc*x + Qg_rc@yintercept
res@tiYAxisString  = f->Qg_tavg@long_name+" ("+f->Qg_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qg,res)

; Evap_rc    = regline(x, Evap)          ; slope
; Evap_arry  = new ( (/2,time_tot/), typeof(Evap))
; Evap_arry(0,:) = Evap
; Evap_arry(1,:) = Evap_rc*x + Evap_rc@yintercept
res@tiYAxisString  = f->Evap_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Evap,res)

; TVeg_rc    = regline(x, TVeg)          ; slope
; TVeg_arry  = new ( (/2,time_tot/), typeof(TVeg))
; TVeg_arry(0,:) = TVeg
; TVeg_arry(1,:) = TVeg_rc*x + TVeg_rc@yintercept
res@tiYAxisString  = f->TVeg_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,TVeg,res)

; ESoil_rc    = regline(x, ESoil)          ; slope
; ESoil_arry  = new ( (/2,time_tot/), typeof(ESoil))
; ESoil_arry(0,:) = ESoil
; ESoil_arry(1,:) = ESoil_rc*x + ESoil_rc@yintercept
res@tiYAxisString  = f->ESoil_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,ESoil,res)

; Qs_rc    = regline(x, Qs)          ; slope
; Qs_arry  = new ( (/2,time_tot/), typeof(Qs))
; Qs_arry(0,:) = Qs
; Qs_arry(1,:) = Qs_rc*x + Qs_rc@yintercept
res@tiYAxisString  = f->Qs_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qs,res)

; Qsb_rc    = regline(x, Qsb)          ; slope
; Qsb_arry  = new ( (/2,time_tot/), typeof(Qsb))
; Qsb_arry(0,:) = Qsb
; Qsb_arry(1,:) = Qsb_rc*x + Qsb_rc@yintercept
res@tiYAxisString  = f->Qsb_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qsb,res)

; SnowDepth_rc    = regline(x, SnowDepth)          ; slope
; SnowDepth_arry  = new ( (/2,time_tot/), typeof(SnowDepth))
; SnowDepth_arry(0,:) = SnowDepth
; SnowDepth_arry(1,:) = SnowDepth_rc*x + SnowDepth_rc@yintercept
res@tiYAxisString  = f->SnowDepth_inst@long_name+" ("+f->SnowDepth_inst@units+")"
plot = gsn_csm_xy(wks,x,SnowDepth,res)

do soil = 0,5
  res@trYMaxF = 50.
  res@trYMinF = 0.
  ; SoilTemp_rc    = regline(x, SoilTemp(:,soil))          ; slope
  ; SoilTemp_arry  = new ( (/2,time_tot/), typeof(SoilTemp))
  ; SoilTemp_arry(0,:) = SoilTemp(:,soil)
  ; SoilTemp_arry(1,:) = SoilTemp_rc*x + SoilTemp_rc@yintercept
  res@tiYAxisString  = "Layer "+(soil+1)+" Soil Temp (C)"
  plot = gsn_csm_xy(wks,x,SoilTemp,res)
  ; delete(SoilTemp_rc)
  ; delete(SoilTemp_arry)
  delete(res@trYMaxF)
  delete(res@trYMinF)
end do
do soil = 0,5
  res@trYMaxF = 0.4
  res@trYMinF = 0.
  ; SoilMoist_rc    = regline(x, SoilMoist(:,soil))          ; slope
  ; SoilMoist_arry  = new ( (/2,time_tot/), typeof(SoilMoist))
  ; SoilMoist_arry(0,:) = SoilMoist(:,soil)
  ; SoilMoist_arry(1,:) = SoilMoist_rc*x + SoilMoist_rc@yintercept
  res@tiYAxisString  = "Layer "+(soil+1)+" Soil Moisture (m3/m3)"
  plot = gsn_csm_xy(wks,x,SoilMois,res)
  ; delete(SoilMoist_rc)
  ; delete(SoilMoist_arry)
  delete(res@trYMaxF)
  delete(res@trYMinF)
end do

res@trYMaxF = 50.
res@trYMinF = 0.
; AvgSurfT_rc    = regline(x, AvgSurfT)          ; slope
; AvgSurfT_arry  = new ( (/2,time_tot/), typeof(AvgSurfT))
; AvgSurfT_arry(0,:) = AvgSurfT
; AvgSurfT_arry(1,:) = AvgSurfT_rc*x + AvgSurfT_rc@yintercept
res@tiYAxisString  = f->AvgSurfT_tavg@long_name+" (C)"
plot = gsn_csm_xy(wks,x,AvgSurfT,res)
delete(res@trYMaxF)
delete(res@trYMinF)

; SoilWet_rc    = regline(x, SoilWet)          ; slope
; SoilWet_arry  = new ( (/2,time_tot/), typeof(SoilWet))
; SoilWet_arry(0,:) = SoilWet
; SoilWet_arry(1,:) = SoilWet_rc*x + SoilWet_rc@yintercept
; res@tiYAxisString = f->SoilWet_inst@long_name+" ("+f->SoilWet_inst@units+")"
; plot = gsn_csm_xy(wks,x,SoilWet_arry,res)

end
