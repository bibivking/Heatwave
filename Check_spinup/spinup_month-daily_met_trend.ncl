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
year_s     = 1970
year_e     = 1999
year_sum   = year_e-year_s+1
pic_name   = "CABLE-HESS_MET_TREND_PFT-"+PFT+"_"+time_scale+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL/"

case_name = (/"outputs-spinup200yr"/)
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
  time_tot = year_sum*365 + 1
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

filename1 = path(0)+"/cable_out_1990_SE_Aus.nc" ; /outputs
f1 = addfile (filename1,"r")
PFT_mask = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
default_val = f1->Qle@_FillValue

delete(filename1)
delete(f1)

Rain     = new((/case_sum,time_tot/),float,default_val)
Tair     = new((/case_sum,time_tot/),float,default_val)
SWdown   = new((/case_sum,time_tot/),float,default_val)
LWdown   = new((/case_sum,time_tot/),float,default_val)
Qair     = new((/case_sum,time_tot/),float,default_val)
Wind     = new((/case_sum,time_tot/),float,default_val)
PSurf    = new((/case_sum,time_tot/),float,default_val)

Rain     = 0.
Tair     = 0.
SWdown   = 0.
LWdown   = 0.
Qair     = 0.
Wind     = 0.
PSurf    = 0.

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

    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc" ; outputs/
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

    SWdown(case_num,i:i+ts-1)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->SWdown(:,{A:B},{C:D}),default_val) ,(/1,2/))
    LWdown(case_num,i:i+ts-1)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->LWdown(:,{A:B},{C:D}),default_val) ,(/1,2/))
    Qair(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qair(:,{A:B},{C:D}),default_val) ,(/1,2/))
    Tair(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Tair(:,{A:B},{C:D}),default_val) ,(/1,2/)) -273.15
    Wind(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Wind(:,{A:B},{C:D}),default_val) ,(/1,2/))
    PSurf(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->PSurf(:,{A:B},{C:D}),default_val) ,(/1,2/))
    Rain(case_num,i:i+ts-1)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Rainf(:,{A:B},{C:D}),default_val) ,(/1,2/))*3600*24

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
Rain_rc    = regline(x, Rain(0,:))          ; slope
Rain_arry  = new ( (/2,time_tot/), typeof(Rain))
Rain_arry(0,:) = Rain(0,:)
Rain_arry(1,:) = Rain_rc*x + Rain_rc@yintercept
res@tiYAxisString  = "Rainfall (mm/d)"
res@tiXAxisString  = "months"
plot = gsn_csm_xy(wks,x,Rain_arry,res)

print("plot Tair")
Tair_rc    = regline(x, Tair(0,:))          ; slope
Tair_arry  = new ( (/2,time_tot/), typeof(Tair))
Tair_arry(0,:) = Tair(0,:)
Tair_arry(1,:) = Tair_rc*x + Tair_rc@yintercept
res@tiYAxisString  = "Tair (C)"
plot = gsn_csm_xy(wks,x,Tair_arry,res)

print("plot Qair")
Qair_rc    = regline(x, Qair(0,:))          ; slope
Qair_arry  = new ( (/2,time_tot/), typeof(Qair))
Qair_arry(0,:) = Qair(0,:)
Qair_arry(1,:) = Qair_rc*x + Qair_rc@yintercept
res@tiYAxisString  = "Qair (kg/kg)"
plot = gsn_csm_xy(wks,x,Qair_arry,res)

print("plot Wind")
Wind_rc    = regline(x, Wind(0,:))          ; slope
Wind_arry  = new ( (/2,time_tot/), typeof(Wind))
Wind_arry(0,:) = Wind(0,:)
Wind_arry(1,:) = Wind_rc*x + Wind_rc@yintercept
res@tiYAxisString  = "Wind (m/s)"
plot = gsn_csm_xy(wks,x,Wind_arry,res)

print("plot SWdown")
SWdown_rc    = regline(x, SWdown(0,:))          ; slope
SWdown_arry  = new ( (/2,time_tot/), typeof(SWdown))
SWdown_arry(0,:) = SWdown(0,:)
SWdown_arry(1,:) = SWdown_rc*x + SWdown_rc@yintercept
res@tiYAxisString  = "SWdown (W/m2)"
plot = gsn_csm_xy(wks,x,SWdown_arry,res)

print("plot LWdown")
LWdown_rc    = regline(x, LWdown(0,:))          ; slope
LWdown_arry  = new ( (/2,time_tot/), typeof(LWdown))
LWdown_arry(0,:) = LWdown(0,:)
LWdown_arry(1,:) = LWdown_rc*x + LWdown_rc@yintercept
res@tiYAxisString  = "LWdown (W/m2)"
plot = gsn_csm_xy(wks,x,LWdown_arry,res)

print("plot PSurf")
PSurf_rc    = regline(x, PSurf(0,:))          ; slope
PSurf_arry  = new ( (/2,time_tot/), typeof(PSurf))
PSurf_arry(0,:) = PSurf(0,:)
PSurf_arry(1,:) = PSurf_rc*x + PSurf_rc@yintercept
res@tiYAxisString  = "PSurf (hPa)"
plot = gsn_csm_xy(wks,x,PSurf_arry,res)

end
