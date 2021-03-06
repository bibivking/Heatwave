;*******************************************
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

path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ERAI_ctl_watmove_new_veg/LIS_output/"
path2 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ERAI_ctl_watmove_new_veg_Gwwb03/LIS_output/"


Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; ____________________________ OPTIONS _______________________________
hydr_param = True
scale      = "AU"
case_name  = "ERAI_watmove"

; region: AU
A = -44.
B = -10.
C = 112.
D = 154.


;______________________________ Data _________________________________
year_s = 1980
year_e = 1980
year_sum = year_e-year_s+1
time_tot = year_sum*365 + 1;0

; Rain     = new((/2,time_tot/),float)
Qle      = new((/2,time_tot/),float)
Qh       = new((/2,time_tot/),float)
Qg       = new((/2,time_tot/),float)
Evap     = new((/2,time_tot/),float)
TVeg     = new((/2,time_tot/),float)
ESoil    = new((/2,time_tot/),float)
Qs       = new((/2,time_tot/),float)
Qsb      = new((/2,time_tot/),float)
Qrec     = new((/2,time_tot/),float)
; AvgSurfT = new((/2,time_tot/),float)

if hydr_param then
  GWwb     = new((/2,time_tot/),float)
  ; Fwsoil   = new((/2,time_tot/),float)
  WatTab   = new((/2,time_tot/),float)
else
  SnowDepth= new((/2,time_tot/),float)
end if

; SoilTemp = new((/2,time_tot,6/),float)
SoilMoist= new((/2,time_tot,6/),float)

; Rain     = 0.
Qle      = 0.
Qh       = 0.
Qg       = 0.
Evap     = 0.
TVeg     = 0.
ESoil    = 0.
Qs       = 0.
Qsb      = 0.
Qrec     = 0.
; AvgSurfT = 0.

if hydr_param then
  GWwb     = 0.
  ; Fwsoil   = 0.
  WatTab   = 0.
else
  SnowDepth= 0.
end if

; SoilTemp = 0.
SoilMoist= 0.

i = 0
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
        filename1 = path2+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
        filename1 = path2+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      f = addfile (filename,"r")
      f1 = addfile (filename1,"r")

      ; ========================= READING =============================

      do day = 1, dom(month-1)
        Qle(0,i)      = avg(f->Qle_tavg(day-1,:,:))
        Qh(0,i)       = avg(f->Qh_tavg(day-1,:,:))
        Qg(0,i)       = avg(f->Qg_tavg(day-1,:,:))
        Evap(0,i)     = avg(f->Evap_tavg(day-1,:,:))*3600*24
        TVeg(0,i)     = avg(f->TVeg_tavg(day-1,:,:))*3600*24
        ESoil(0,i)    = avg(f->ESoil_tavg(day-1,:,:))*3600*24
        Qs(0,i)       = avg(f->Qs_tavg(day-1,:,:))*3600*24
        Qsb(0,i)      = avg(f->Qsb_tavg(day-1,:,:))*3600*24
        Qrec(0,i)     = avg(f->Qrec_tavg(day-1,:,:))*3600*24

        Qle(1,i)      = avg(f1->Qle_tavg(day-1,:,:))
        Qh(1,i)       = avg(f1->Qh_tavg(day-1,:,:))
        Qg(1,i)       = avg(f1->Qg_tavg(day-1,:,:))
        Evap(1,i)     = avg(f1->Evap_tavg(day-1,:,:))*3600*24
        TVeg(1,i)     = avg(f1->TVeg_tavg(day-1,:,:))*3600*24
        ESoil(1,i)    = avg(f1->ESoil_tavg(day-1,:,:))*3600*24
        Qs(1,i)       = avg(f1->Qs_tavg(day-1,:,:))*3600*24
        Qsb(1,i)      = avg(f1->Qsb_tavg(day-1,:,:))*3600*24
        Qrec(1,i)     = avg(f1->Qrec_tavg(day-1,:,:))*3600*24

        ; AvgSurfT(0,i) = avg(f->AvgSurfT_tavg(day-1,:,:))-273.15
        ; AvgSurfT(1,i) = avg(f1->AvgSurfT_tavg(day-1,:,:))-273.15

        if hydr_param then
          GWwb(0,i)    = avg(f->GWwb_inst(day-1,:,:))*25
          ; Fwsoil(0,i)  = avg(f->fwsoil_inst(day-1,:,:))
          WatTab(0,i)  = avg(f->WaterTableD_tavg(day-1,:,:))/1000.
          ; Rain(0,i)    = avg(f->TotalPrecip_tavg(day-1,:,:))*3600*24
          GWwb(1,i)    = avg(f1->GWwb_inst(day-1,:,:))*25
          ; Fwsoil(1,i)  = avg(f1->fwsoil_inst(day-1,:,:))
          WatTab(1,i)  = avg(f1->WaterTableD_tavg(day-1,:,:))/1000.
          ; Rain(1,i)    = avg(f1->TotalPrecip_tavg(day-1,:,:))*3600*24
        else
          ; Rain(0,i)     = avg(f->Rainf_f_inst(day-1,:,:))*3600*24
          SnowDepth(0,i)= avg(f->SnowDepth_inst(day-1,:,:))
          ; Rain(1,i)     = avg(f1->Rainf_f_inst(day-1,:,:))*3600*24
          SnowDepth(1,i)= avg(f1->SnowDepth_inst(day-1,:,:))
        end if

        do soil = 0,5
          ; SoilTemp(0,i,soil) = avg(f->SoilTemp_inst(day-1,soil,:,:))-273.15
          SoilMoist(0,i,soil)= avg(f->SoilMoist_inst(day-1,soil,:,:))
          ; SoilTemp(1,i,soil) = avg(f1->SoilTemp_inst(day-1,soil,:,:))-273.15
          SoilMoist(1,i,soil)= avg(f1->SoilMoist_inst(day-1,soil,:,:))
        end do
        i = i + 1
      end do ; day
      ; delete(filename)
      ; delete(f)
    end do ; month
    ;delete(f1)
end do ; year
print("prepare to plot")


;_____________________________ PLOT ___________________________________
if hydr_param then
  pic = "./plots/Spinup_LIS_vs_CABLE-HESS-"+case_name+"_"+scale
else
  pic = "./plots/Spinup_CABLE-HESS_GSWP3"
end if

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
res@pmLegendDisplayMode = "Always"
res@xyLineColors       = (/"blue","red","green","orange","red","red","red"/)          ; change line color
res@xyExplicitLegendLabels = (/"GWwb=0.45","GWwb=0.3","GW_1-60","NoAQ"/)
;res@trYMinF = 0


x             = ispan(1,time_tot,1)
res@tmXBMode  = "Explicit"
res@tmXBValues= ispan(1,time_tot+1,365*5)
res@tmXBLabels= tostring(ispan(year_s,year_e+1,5))



filename = path1+"/LIS.CABLE.1999090100.d01.nc"
f = addfile(filename,"r")
;
; print("plot rainfall")
; ; Rain_rc    = regline(x, Rain)          ; slope
; ; Rain_arry  = new ( (/2,time_tot/), typeof(Rain))
; ; Rain_arry(0,:) = Rain
; ; Rain_arry(1,:) = Rain_rc*x + Rain_rc@yintercept
; if hydr_param then
;   res@tiYAxisString  = f->TotalPrecip_tavg@long_name+" (mm/d)"
; else
;   res@tiYAxisString  = f->Rainf_f_inst@long_name+" (mm/d)"
; end if
; res@tiXAxisString  = "days"
; plot = gsn_csm_xy(wks,x,Rain,res)
;
print("plot latent heat")
; Qle_rc    = regline(x, Qle)          ; slope
; Qle_arry  = new ( (/2,time_tot/), typeof(Qle))
; Qle_arry(0,:) = Qle
; Qle_arry(1,:) = Qle_rc*x + Qle_rc@yintercept
res@tiYAxisString  = f->Qle_tavg@long_name+" ("+f->Qle_tavg@units+")"
res@tiXAxisString  = "days"
plot = gsn_csm_xy(wks,x,Qle,res)

print("plot sensible heat")
; Qh_rc    = regline(x, Qh)          ; slope
; Qh_arry  = new ( (/2,time_tot/), typeof(Qh))
; Qh_arry(0,:) = Qh
; Qh_arry(1,:) = Qh_rc*x + Qh_rc@yintercept
res@tiYAxisString  = f->Qh_tavg@long_name+" ("+f->Qh_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qh,res)

print("plot ground heat")
; Qg_rc    = regline(x, Qg)          ; slope
; Qg_arry  = new ( (/2,time_tot/), typeof(Qg))
; Qg_arry(0,:) = Qg
; Qg_arry(1,:) = Qg_rc*x + Qg_rc@yintercept
res@tiYAxisString  = f->Qg_tavg@long_name+" ("+f->Qg_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qg,res)

print("plot ET")
; Evap_rc    = regline(x, Evap)          ; slope
; Evap_arry  = new ( (/2,time_tot/), typeof(Evap))
; Evap_arry(0,:) = Evap
; Evap_arry(1,:) = Evap_rc*x + Evap_rc@yintercept
res@tiYAxisString  = f->Evap_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Evap,res)

print("plot trans")
; TVeg_rc    = regline(x, TVeg)          ; slope
; TVeg_arry  = new ( (/2,time_tot/), typeof(TVeg))
; TVeg_arry(0,:) = TVeg
; TVeg_arry(1,:) = TVeg_rc*x + TVeg_rc@yintercept
res@tiYAxisString  = f->TVeg_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,TVeg,res)

print("plot esoil")
; ESoil_rc    = regline(x, ESoil)          ; slope
; ESoil_arry  = new ( (/2,time_tot/), typeof(ESoil))
; ESoil_arry(0,:) = ESoil
; ESoil_arry(1,:) = ESoil_rc*x + ESoil_rc@yintercept
res@tiYAxisString  = f->ESoil_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,ESoil,res)

print("plot surface runoff")
; Qs_rc    = regline(x, Qs)          ; slope
; Qs_arry  = new ( (/2,time_tot/), typeof(Qs))
; Qs_arry(0,:) = Qs
; Qs_arry(1,:) = Qs_rc*x + Qs_rc@yintercept
res@tiYAxisString  = f->Qs_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qs,res)

print("plot subsurface runoff")
; Qsb_rc    = regline(x, Qsb)          ; slope
; Qsb_arry  = new ( (/2,time_tot/), typeof(Qsb))
; Qsb_arry(0,:) = Qsb
; Qsb_arry(1,:) = Qsb_rc*x + Qsb_rc@yintercept
res@tiYAxisString  = f->Qsb_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qsb,res)

print("plot aquifer recharge")
; Qrec_rc    = regline(x, Qrec)          ; slope
; Qrec_arry  = new ( (/2,time_tot/), typeof(Qrec))
; Qrec_arry(0,:) = Qrec
; Qrec_arry(1,:) = Qrec_rc*x + Qrec_rc@yintercept
res@tiYAxisString  = f->Qrec_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qrec,res)


if hydr_param then
  print("plot aquifer moisture")
  ; GWwb_rc    = regline(x, GWwb)          ; slope
  ; GWwb_arry  = new ( (/2,time_tot/), typeof(GWwb))
  ; GWwb_arry(0,:) = GWwb
  ; GWwb_arry(1,:) = GWwb_rc*x + GWwb_rc@yintercept
  res@trYMaxF = 0.5
  res@trYMinF = 0.
  res@tiYAxisString  = "Aquifer water storage (m)"
  plot = gsn_csm_xy(wks,x,GWwb,res)
  delete(res@trYMaxF)
  delete(res@trYMinF)
  ; print("plot fwsoil")
  ; ; fwsoil_rc    = regline(x, Fwsoil)          ; slope
  ; ; fwsoil_arry  = new ( (/2,time_tot/), typeof(Fwsoil))
  ; ; fwsoil_arry(0,:) = Fwsoil
  ; ; fwsoil_arry(1,:) = fwsoil_rc*x + fwsoil_rc@yintercept
  ; res@tiYAxisString  = f->fwsoil_inst@long_name+" ("+f->fwsoil_inst@units+")"
  ; plot = gsn_csm_xy(wks,x,Fwsoil,res)
  print("plot water table")
  ; WatTab_rc    = regline(x, WatTab)          ; slope
  ; WatTab_arry  = new ( (/2,time_tot/), typeof(WatTab))
  ; WatTab_arry(0,:) = WatTab
  ; WatTab_arry(1,:) = WatTab_rc*x + WatTab_rc@yintercept
  res@trYMaxF = 12.
  res@trYMinF = 0.
  res@tiYAxisString  = f->WaterTableD_tavg@long_name+" ("+f->WaterTableD_tavg@units+")"
  plot = gsn_csm_xy(wks,x,WatTab,res)
  delete(res@trYMaxF)
  delete(res@trYMinF)
else
  print("plot snow depth")
  ; SnowDepth_rc    = regline(x, SnowDepth)          ; slope
  ; SnowDepth_arry  = new ( (/2,time_tot/), typeof(SnowDepth))
  ; SnowDepth_arry(0,:) = SnowDepth
  ; SnowDepth_arry(1,:) = SnowDepth_rc*x + SnowDepth_rc@yintercept
  res@tiYAxisString  = f->SnowDepth_inst@long_name+" ("+f->SnowDepth_inst@units+")"
  plot = gsn_csm_xy(wks,x,SnowDepth,res)
end if

;
; do soil = 0,5
;   print("plot soil temp")
;   res@trYMaxF = 50.
;   res@trYMinF = 0.
;   ; SoilTemp_rc    = regline(x, SoilTemp(:,soil))          ; slope
;   ; SoilTemp_arry  = new ( (/2,time_tot/), typeof(SoilTemp))
;   ; SoilTemp_arry(0,:) = SoilTemp(:,soil)
;   ; SoilTemp_arry(1,:) = SoilTemp_rc*x + SoilTemp_rc@yintercept
;   res@tiYAxisString  = "Layer "+(soil+1)+" Soil Temp (C)"
;   plot = gsn_csm_xy(wks,x,SoilTemp(:,:,soil),res)
;   ; delete(SoilTemp_rc)
;   ; delete(SoilTemp_arry)
;   delete(res@trYMaxF)
;   delete(res@trYMinF)
; end do
do soil = 0,5
  print("plot soil moisture")
  res@trYMaxF = 0.4
  res@trYMinF = 0.
  ; SoilMoist_rc    = regline(x, SoilMoist(:,soil))          ; slope
  ; SoilMoist_arry  = new ( (/2,time_tot/), typeof(SoilMoist))
  ; SoilMoist_arry(0,:) = SoilMoist(:,soil)
  ; SoilMoist_arry(1,:) = SoilMoist_rc*x + SoilMoist_rc@yintercept
  res@tiYAxisString  = "Layer "+(soil+1)+" Soil Moisture (m3/m3)"
  plot = gsn_csm_xy(wks,x,SoilMoist(:,:,soil),res)
  ; delete(SoilMoist_rc)
  ; delete(SoilMoist_arry)
  delete(res@trYMaxF)
  delete(res@trYMinF)
end do
;
; print("plot surface temp")
; res@trYMaxF = 50.
; res@trYMinF = 0.
; ; AvgSurfT_rc    = regline(x, AvgSurfT)          ; slope
; ; AvgSurfT_arry  = new ( (/2,time_tot/), typeof(AvgSurfT))
; ; AvgSurfT_arry(0,:) = AvgSurfT
; ; AvgSurfT_arry(1,:) = AvgSurfT_rc*x + AvgSurfT_rc@yintercept
; res@tiYAxisString  = f->AvgSurfT_tavg@long_name+" (C)"
; plot = gsn_csm_xy(wks,x,AvgSurfT,res)
; delete(res@trYMaxF)
; delete(res@trYMinF)

; SoilWet_rc    = regline(x, SoilWet)          ; slope
; SoilWet_arry  = new ( (/2,time_tot/), typeof(SoilWet))
; SoilWet_arry(0,:) = SoilWet
; SoilWet_arry(1,:) = SoilWet_rc*x + SoilWet_rc@yintercept
; res@tiYAxisString = f->SoilWet_inst@long_name+" ("+f->SoilWet_inst@units+")"
; plot = gsn_csm_xy(wks,x,SoilWet_arry,res)

end
