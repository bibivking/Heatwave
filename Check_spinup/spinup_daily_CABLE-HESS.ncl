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

path1 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3/outputs/GSWP/spinup1"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3/outputs/GSWP/spinup2"
path3 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3/outputs/GSWP/spinup3"

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; ____________________________ OPTIONS _______________________________
hydr_param = True
scale      = "SE-AU"
case_name  = "CABLE-HESS_spinup1-3"

; region: AU
if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -27.
  C = 140.
  D = 155.
end if
end if

;______________________________ Data _________________________________
year_s = 1979
year_e = 2010
year_sum = year_e-year_s+1
time_tot = year_sum*365 + 7;8

; Rain     = new((/3,time_tot/),float)
Qle      = new((/3,time_tot/),float)
Qh       = new((/3,time_tot/),float)
Qg       = new((/3,time_tot/),float)
Evap     = new((/3,time_tot/),float)
TVeg     = new((/3,time_tot/),float)
ESoil    = new((/3,time_tot/),float)
Qs       = new((/3,time_tot/),float)
Qsb      = new((/3,time_tot/),float)
Qrec     = new((/3,time_tot/),float)
; AvgSurfT = new((/3,time_tot/),float)

if hydr_param then
  GWwb     = new((/3,time_tot/),float)
  ; Fwsoil   = new((/3,time_tot/),float)
  WatTab   = new((/3,time_tot/),float)
else
  SnowDepth= new((/3,time_tot/),float)
end if

; SoilTemp = new((/3,time_tot,6/),float)
SoilMoist= new((/3,time_tot,6/),float)

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

    ; read CABLE (HESS)
    filename1 = path1+"/cable_out_"+year+".nc"
    f1 = addfile (filename1,"r")
    filename2 = path2+"/cable_out_"+year+".nc"
    f2 = addfile (filename2,"r")
    filename3 = path3+"/cable_out_"+year+".nc"
    f3 = addfile (filename3,"r")

    day_acu = 1

    ; read LIS-CABLE
    do month = 1, 12

      do day = 1, dom(month-1)

        ; --------------- read spinup1 ------------------
        Qle(0,i)      = avg(f1->Qle(day_acu-1,{A:B},{C:D}))
        Qh(0,i)       = avg(f1->Qh(day_acu-1,{A:B},{C:D}))
        Qg(0,i)       = avg(f1->Qg(day_acu-1,{A:B},{C:D}))
        Evap(0,i)     = avg(f1->Evap(day_acu-1,{A:B},{C:D}))*3600*24
        TVeg(0,i)     = avg(f1->TVeg(day_acu-1,{A:B},{C:D}))*3600*24
        ESoil(0,i)    = avg(f1->ESoil(day_acu-1,{A:B},{C:D}))*3600*24
        Qs(0,i)       = avg(f1->Qs(day_acu-1,{A:B},{C:D}))*3600*24
        Qsb(0,i)      = avg(f1->Qsb(day_acu-1,{A:B},{C:D}))*3600*24
        Qrec(0,i)     = avg(f1->Qrecharge(day_acu-1,{A:B},{C:D}))*3600*24

        if hydr_param then
          GWwb(0,i)    = avg(f1->GWMoist(day_acu-1,{A:B},{C:D}))
          WatTab(0,i)  = avg(f1->WatTable(day_acu-1,{A:B},{C:D}))
        else
          SnowDepth(0,i)= avg(f1->SnowDepth(day_acu-1,{A:B},{C:D}))
        end if

        do soil = 0,5
        ;   SoilTemp(0,i,soil) = avg(f1->SoilTemp(day_acu-1,soil,{A:B},{C:D}))-273.15
          SoilMoist(0,i,soil)= avg(f1->SoilMoist(day_acu-1,soil,{A:B},{C:D}))
        end do

        ;
        ; ; --------------- read CABLE(HESS) ------------------
        ; Qle(1,i)      = avg(f2->Qle(day_acu-1,{A:B},{C:D}))
        ; Qh(1,i)       = avg(f2->Qh(day_acu-1,{A:B},{C:D}))
        ; Qg(1,i)       = avg(f2->Qg(day_acu-1,{A:B},{C:D}))
        ; Evap(1,i)     = avg(f2->Evap(day_acu-1,{A:B},{C:D}))*3600*24
        ; TVeg(1,i)     = avg(f2->TVeg(day_acu-1,{A:B},{C:D}))*3600*24
        ; ESoil(1,i)    = avg(f2->ESoil(day_acu-1,{A:B},{C:D}))*3600*24
        ; Qs(1,i)       = avg(f2->Qs(day_acu-1,{A:B},{C:D}))*3600*24
        ; Qsb(1,i)      = avg(f2->Qsb(day_acu-1,{A:B},{C:D}))*3600*24
        ; Qrec(1,i)     = avg(f2->Qrecharge(day_acu-1,{A:B},{C:D}))*3600*24
        ;
        ; if hydr_param then
        ;   GWwb(1,i)    = avg(f2->GWMoist(day_acu-1,{A:B},{C:D}))
        ;   WatTab(1,i)  = avg(f2->WatTable(day_acu-1,{A:B},{C:D}))
        ; else
        ;   SnowDepth(1,i)= avg(f2->SnowDepth(day_acu-1,{A:B},{C:D}))
        ; end if
        ;
        ; do soil = 0,5
        ; ;   SoilTemp(1,i,soil) = avg(f2->SoilTemp(day_acu-1,soil,{A:B},{C:D}))-273.15
        ;   SoilMoist(1,i,soil)= avg(f2->SoilMoist(day_acu-1,soil,{A:B},{C:D}))
        ; end do
        ;
        ; ; --------------- read CABLE(HESS) ------------------
        ; Qle(2,i)      = avg(f3->Qle(day_acu-1,{A:B},{C:D}))
        ; Qh(2,i)       = avg(f3->Qh(day_acu-1,{A:B},{C:D}))
        ; Qg(2,i)       = avg(f3->Qg(day_acu-1,{A:B},{C:D}))
        ; Evap(2,i)     = avg(f3->Evap(day_acu-1,{A:B},{C:D}))*3600*24
        ; TVeg(2,i)     = avg(f3->TVeg(day_acu-1,{A:B},{C:D}))*3600*24
        ; ESoil(2,i)    = avg(f3->ESoil(day_acu-1,{A:B},{C:D}))*3600*24
        ; Qs(2,i)       = avg(f3->Qs(day_acu-1,{A:B},{C:D}))*3600*24
        ; Qsb(2,i)      = avg(f3->Qsb(day_acu-1,{A:B},{C:D}))*3600*24
        ; Qrec(2,i)     = avg(f3->Qrecharge(day_acu-1,{A:B},{C:D}))*3600*24
        ;
        ; if hydr_param then
        ;   GWwb(2,i)    = avg(f3->GWMoist(day_acu-1,{A:B},{C:D}))
        ;   WatTab(2,i)  = avg(f3->WatTable(day_acu-1,{A:B},{C:D}))
        ; else
        ;   SnowDepth(2,i)= avg(f3->SnowDepth(day_acu-1,{A:B},{C:D}))
        ; end if
        ;
        ; do soil = 0,5
        ; ;   SoilTemp(2,i,soil) = avg(f3->SoilTemp(day_acu-1,soil,{A:B},{C:D}))-273.15
        ;   SoilMoist(2,i,soil)= avg(f3->SoilMoist(day_acu-1,soil,{A:B},{C:D}))
        ; end do

        day_acu = day_acu + 1
        i = i + 1
      end do ; day
    end do ; month
end do ; year
print("prepare to plot")


;_____________________________ PLOT ___________________________________
if hydr_param then
  pic = "./plots/Spinup_"+case_name+"_"+scale
else
  pic = "./plots/Spinup_"+case_name+"_"+scale
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
res@xyExplicitLegendLabels = (/"spinup1","spinup2","spinup3","NoAQ"/)
;res@trYMinF = 0


x             = ispan(1,time_tot,1)
res@tmXBMode  = "Explicit"
res@tmXBValues= ispan(1,time_tot+1,365*5)
res@tmXBLabels= tostring(ispan(year_s,year_e+1,5))

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
res@tiYAxisString  = "Latent Heat Flux (W/m2)"
res@tiXAxisString  = "days"
plot = gsn_csm_xy(wks,x,Qle,res)

print("plot sensible heat")
; Qh_rc    = regline(x, Qh)          ; slope
; Qh_arry  = new ( (/2,time_tot/), typeof(Qh))
; Qh_arry(0,:) = Qh
; Qh_arry(1,:) = Qh_rc*x + Qh_rc@yintercept
res@tiYAxisString  = "Sensible Heat Flux (W/m2)"
plot = gsn_csm_xy(wks,x,Qh,res)

print("plot ground heat")
; Qg_rc    = regline(x, Qg)          ; slope
; Qg_arry  = new ( (/2,time_tot/), typeof(Qg))
; Qg_arry(0,:) = Qg
; Qg_arry(1,:) = Qg_rc*x + Qg_rc@yintercept
res@tiYAxisString  = "Ground Heat Flux (W/m2)"
plot = gsn_csm_xy(wks,x,Qg,res)

print("plot ET")
; Evap_rc    = regline(x, Evap)          ; slope
; Evap_arry  = new ( (/2,time_tot/), typeof(Evap))
; Evap_arry(0,:) = Evap
; Evap_arry(1,:) = Evap_rc*x + Evap_rc@yintercept
res@tiYAxisString  = "Total ET (mm/d)"
plot = gsn_csm_xy(wks,x,Evap,res)

print("plot trans")
; TVeg_rc    = regline(x, TVeg)          ; slope
; TVeg_arry  = new ( (/2,time_tot/), typeof(TVeg))
; TVeg_arry(0,:) = TVeg
; TVeg_arry(1,:) = TVeg_rc*x + TVeg_rc@yintercept
res@tiYAxisString  = "Transpiration (mm/d)"
plot = gsn_csm_xy(wks,x,TVeg,res)

print("plot esoil")
; ESoil_rc    = regline(x, ESoil)          ; slope
; ESoil_arry  = new ( (/2,time_tot/), typeof(ESoil))
; ESoil_arry(0,:) = ESoil
; ESoil_arry(1,:) = ESoil_rc*x + ESoil_rc@yintercept
res@tiYAxisString  = "Soil Evaporation (mm/d)"
plot = gsn_csm_xy(wks,x,ESoil,res)

print("plot surface runoff")
; Qs_rc    = regline(x, Qs)          ; slope
; Qs_arry  = new ( (/2,time_tot/), typeof(Qs))
; Qs_arry(0,:) = Qs
; Qs_arry(1,:) = Qs_rc*x + Qs_rc@yintercept
res@tiYAxisString  = "Surface Runoff (mm/d)"
plot = gsn_csm_xy(wks,x,Qs,res)

print("plot subsurface runoff")
; Qsb_rc    = regline(x, Qsb)          ; slope
; Qsb_arry  = new ( (/2,time_tot/), typeof(Qsb))
; Qsb_arry(0,:) = Qsb
; Qsb_arry(1,:) = Qsb_rc*x + Qsb_rc@yintercept
res@tiYAxisString  = "Subsurface Runoff (mm/d)"
plot = gsn_csm_xy(wks,x,Qsb,res)


if hydr_param then
  print("plot aquifer moisture")
  ; GWwb_rc    = regline(x, GWwb)          ; slope
  ; GWwb_arry  = new ( (/2,time_tot/), typeof(GWwb))
  ; GWwb_arry(0,:) = GWwb
  ; GWwb_arry(1,:) = GWwb_rc*x + GWwb_rc@yintercept
  res@trYMaxF = 0.5
  res@trYMinF = 0.
  res@tiYAxisString  = "Aquifer Moisture (m3/m3)"
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
  res@trYMaxF = 20.
  res@trYMinF = 0.
  res@tiYAxisString  = "Water Table Depth (m)"
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
