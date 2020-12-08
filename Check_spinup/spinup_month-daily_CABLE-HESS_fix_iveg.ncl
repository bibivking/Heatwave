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

; ____________________________ OPTIONS _______________________________
hydr_param = True
time_scale = "monthly" ;"monthly" ; "daily"
PFT        = "all" ; "tree"
scale      = "SE-AU"
year_s     = 1971
year_e     = 1971
year_sum   = year_e-year_s+1
pic_name   = "CABLE-HESS_Fix-iveg_PFT-"+PFT+"_"+time_scale+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"

case_name = (/"SPINUP_30yrx2","SPINUP_30yrx2_fix_iveg"/)
Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers


if PFT .eq. "all" then
  iveg_s = 0
  iveg_e = 20
else if PFT .eq. "tree" then
  iveg_s = 1
  iveg_e = 4
end if
end if

if time_scale .eq. "monthly" then
  print("monthly")
  time_tot = year_sum*12
else if time_scale .eq. "daily" then
  print("daily")
  time_tot = year_sum*365 
  print(time_tot)
end if
end if

; region: AU
if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -28.
  C = 140.
  D = 154.
end if
end if

;______________________________ Data _________________________________
case_sum  = dimsizes(case_name)
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = path_mum+case_name(case_num)
end do

; ___________________________ TREE Mask ______________________________

filename1 = path(0)+"/outputs-second/cable_out_1971_SE_Aus.nc" ; /outputs
f1 = addfile (filename1,"r")
PFT_mask = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
default_val = f1->Qle@_FillValue

delete(filename1)
delete(f1)

Rain     = new((/case_sum,time_tot/),float,default_val)
Qle      = new((/case_sum,time_tot/),float,default_val)
Qh       = new((/case_sum,time_tot/),float,default_val)
Qg       = new((/case_sum,time_tot/),float,default_val)
Evap     = new((/case_sum,time_tot/),float,default_val)
TVeg     = new((/case_sum,time_tot/),float,default_val)
ESoil    = new((/case_sum,time_tot/),float,default_val)
Qs       = new((/case_sum,time_tot/),float,default_val)
Qsb      = new((/case_sum,time_tot/),float,default_val)
Qrec     = new((/case_sum,time_tot/),float,default_val)
VegT     = new((/case_sum,time_tot/),float,default_val)
CanT     = new((/case_sum,time_tot/),float,default_val)

if hydr_param then
  GWwb     = new((/case_sum,time_tot/),float,default_val)
  Fwsoil   = new((/case_sum,time_tot/),float,default_val)
  WatTab   = new((/case_sum,time_tot/),float,default_val)
else
  SnowDepth= new((/case_sum,time_tot/),float,default_val)
end if

SoilTemp = new((/case_sum,time_tot,6/),float,default_val)
SoilMoist= new((/case_sum,time_tot,6/),float,default_val)

Rain     = 0.
Qle      = 0.
Qh       = 0.
Qg       = 0.
Evap     = 0.
TVeg     = 0.
ESoil    = 0.
Qs       = 0.
Qsb      = 0.
Qrec     = 0.
VegT     = 0.
CanT     = 0.

if hydr_param then
  GWwb     = 0.
  Fwsoil   = 0.
  WatTab   = 0.
else
  SnowDepth= 0.
end if

SoilTemp = 0.
SoilMoist= 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  i = 0
  do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if
    if case_num .eq. 0 then
       filename1 = path(case_num)+"/outputs-second/cable_out_"+year+"_SE_Aus.nc" ; outputs/
    else 
       filename1 = path(case_num)+"/outputs/cable_out_"+year+"_SE_Aus.nc"
    end if
    f1 = addfile (filename1,"r")

    if time_scale .eq. "monthly" then
      ts = 12
    else if time_scale .eq. "daily" then
      ts = sum(dom)
    end if
    end if

    PFT_mask_multi = new((/ts,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
    do z = 0, ts-1
      PFT_mask_multi(z,:,:) = PFT_mask
    end do

    Qle(case_num,i:i+ts-1)      = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qle(:,{A:B},{C:D}),default_val) ,(/1,2/))
    Qh(case_num,i:i+ts-1)       = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qh(:,{A:B},{C:D}),default_val) ,(/1,2/))
    Qg(case_num,i:i+ts-1)       = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qg(:,{A:B},{C:D}),default_val) ,(/1,2/))
    Rain(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Rainf(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    Evap(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Evap(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    TVeg(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->TVeg(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    ESoil(case_num,i:i+ts-1)    = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->ESoil(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    Qs(case_num,i:i+ts-1)       = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qs(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    Qsb(case_num,i:i+ts-1)      = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qsb(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    Qrec(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qrecharge(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24
    VegT(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->VegT(:,{A:B},{C:D}),default_val) ,(/1,2/))-273.15
    CanT(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->CanT(:,{A:B},{C:D}),default_val) ,(/1,2/))-273.15

    if hydr_param then
      GWwb(case_num,i:i+ts-1)    = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1,f1->GWMoist(:,{A:B},{C:D}),default_val) ,(/1,2/))
      Fwsoil(case_num,i:i+ts-1)  = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1,f1->Fwsoil(:,{A:B},{C:D}),default_val) ,(/1,2/))
      WatTab(case_num,i:i+ts-1)  = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1,f1->WatTable(:,{A:B},{C:D}),default_val) ,(/1,2/))
    else
      SnowDepth(case_num,i:i+ts-1)  = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1,f1->SnowDepth(:,{A:B},{C:D}),default_val) ,(/1,2/))
    end if
    do soil = 0,5
       SoilTemp(case_num,i:i+ts-1,soil) = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->SoilTemp(:,soil,{A:B},{C:D}),default_val), (/1,2/))-273.15
       SoilMoist(case_num,i:i+ts-1,soil)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->SoilMoist(:,soil,{A:B},{C:D}),default_val), (/1,2/))
    end do
    i = i + ts
    delete(PFT_mask_multi)
  end do ; year
end do ; case_num
print("prepare to plot")


;_____________________________ PLOT ___________________________________
pic = "./plots/Spinup_"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��
res                    = True
res@gsnMaximize        = True
res@vpHeightF          = 0.5
res@vpWidthF           = 0.90
res@vpXF               = 0.08
res@vpYF               = 0.7
res@pmLegendWidthF     = 0.10                     ; Resize legend width
res@pmLegendHeightF    = 0.10                     ; and height
res@tiXAxisFontHeightF = 0.02
res@tiYAxisFontHeightF = 0.02
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
res@gsnPaperOrientation = "landscape";"portrait";"landscape" ;"auto"

res@xyLineThicknesses   = (/  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0/)          ; make second line thicker
res@pmLegendDisplayMode = "Always"
res@xyDashPatterns      = (/0, 0, 0, 0, 0, 0, 0, 0/)
res@xyLineColors        = (/"black","red","darkgreen","purple","blue","chartreuse","orange","cyan"/)          ; change line color
res@xyExplicitLegendLabels = case_name
;res@trYMinF = 0

x             = ispan(1,time_tot,1)
res@tmXBMode  = "Explicit"
if time_scale .eq. "monthly" then
  res@tmXBValues= ispan(1,time_tot,12)
  res@tiXAxisString  = "months"
else if time_scale .eq. "daily" then
  res@tmXBValues= ispan(1,time_tot,366)
  res@tiXAxisString  = "days"
end if
end if
res@tmXBLabels= tostring(ispan(year_s,year_e,1))


print("plot rainfall")
; Rain_rc    = regline(x, Rain)          ; slope
; Rain_arry  = new ( (/2,time_tot/), typeof(Rain))
; Rain_arry(0,:) = Rain
; Rain_arry(1,:) = Rain_rc*x + Rain_rc@yintercept
res@tiYAxisString  = "Rainfall (mm/d)"
res@tiXAxisString  = "months"
plot = gsn_csm_xy(wks,x,Rain(0,:),res)
;
print("plot latent heat")
; Qle_rc    = regline(x, Qle)          ; slope
; Qle_arry  = new ( (/2,time_tot/), typeof(Qle))
; Qle_arry(0,:) = Qle
; Qle_arry(1,:) = Qle_rc*x + Qle_rc@yintercept
res@tiYAxisString  = "Latent Heat Flux (W/m2)"
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

print("plot aquifer recharge")
; Qrec_rc    = regline(x, Qrec)          ; slope
; Qrec_arry  = new ( (/2,time_tot/), typeof(Qrec))
; Qrec_arry(0,:) = Qrec
; Qrec_arry(1,:) = Qrec_rc*x + Qrec_rc@yintercept
res@tiYAxisString  = "Aquifer Recharge (mm/d)"
plot = gsn_csm_xy(wks,x,Qrec,res)

print("plot vegetation temperature")
; VegT_rc    = regline(x, VegT)          ; slope
; VegT_arry  = new ( (/2,time_tot/), typeof(VegT))
; VegT_arry(0,:) = VegT
; VegT_arry(1,:) = VegT_rc*x + VegT_rc@yintercept
res@tiYAxisString  = "Vegetation Temperature (C)"
plot = gsn_csm_xy(wks,x,VegT,res)

print("plot within-canopy temperature")
; CanT_rc    = regline(x, CanT)          ; slope
; CanT_arry  = new ( (/2,time_tot/), typeof(CanT))
; CanT_arry(0,:) = CanT
; CanT_arry(1,:) = CanT_rc*x + CanT_rc@yintercept
res@tiYAxisString  = "Within-canopy Temperature (C)"
plot = gsn_csm_xy(wks,x,CanT,res)

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
  res@trYMaxF = 1.
  res@trYMinF = 0.
  res@tiYAxisString  = " Fwsoil (-)"
  plot = gsn_csm_xy(wks,x,Fwsoil,res)
  delete(res@trYMaxF)
  delete(res@trYMinF)
  print("plot water table")
  ; WatTab_rc    = regline(x, WatTab)          ; slope
  ; WatTab_arry  = new ( (/2,time_tot/), typeof(WatTab))
  ; WatTab_arry(0,:) = WatTab
  ; WatTab_arry(1,:) = WatTab_rc*x + WatTab_rc@yintercept
  ; res@trYMaxF = 100.
  ; res@trYMinF = 0.
  res@tiYAxisString  = "Water Table Depth (m)"
  plot = gsn_csm_xy(wks,x,WatTab,res)
  ; delete(res@trYMaxF)
  ; delete(res@trYMinF)
else
  print("plot snow depth")
  ; SnowDepth_rc    = regline(x, SnowDepth)          ; slope
  ; SnowDepth_arry  = new ( (/2,time_tot/), typeof(SnowDepth))
  ; SnowDepth_arry(0,:) = SnowDepth
  ; SnowDepth_arry(1,:) = SnowDepth_rc*x + SnowDepth_rc@yintercept
  res@tiYAxisString  = f->SnowDepth_inst@long_name+" ("+f->SnowDepth_inst@units+")"
  plot = gsn_csm_xy(wks,x,SnowDepth,res)
end if


do soil = 0,5
  print("plot soil temp")
  res@trYMaxF = 50.
  res@trYMinF = 0.
  ; SoilTemp_rc    = regline(x, SoilTemp(:,soil))          ; slope
  ; SoilTemp_arry  = new ( (/2,time_tot/), typeof(SoilTemp))
  ; SoilTemp_arry(0,:) = SoilTemp(:,soil)
  ; SoilTemp_arry(1,:) = SoilTemp_rc*x + SoilTemp_rc@yintercept
  res@tiYAxisString  = "Layer "+(soil+1)+" Soil Temp (C)"
  plot = gsn_csm_xy(wks,x,SoilTemp(:,:,soil),res)
  ; delete(SoilTemp_rc)
  ; delete(SoilTemp_arry)
  delete(res@trYMaxF)
  delete(res@trYMinF)
end do
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
end
