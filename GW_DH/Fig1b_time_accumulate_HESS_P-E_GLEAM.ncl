
;MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
; Author  : Mengyuan Mu
; Version : 1.0 (22.03.2021)"
; Email   : mu.mengyuan815@gmail.com
;WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; ================================= Predefine ==================================
pic_name  = "P-E"
case_name = (/"CTL-30x3+30yr/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr"/)

PFT       = "all" ; "tree" ;"nontree"
scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
GLEAM_version = "3.5a";"3.3a"
case_sum  = dimsizes(case_name)

year_s   = 2000
year_e   = 2019
year_sum = year_e-year_s+1

time_tot = 0
do year = year_s,year_e
  if mod(year,4) .eq. 0 then
    time_tot = time_tot + 366
  else
    time_tot = time_tot + 365
  end if
end do


; set PFT
iveg_tree_s    = 1
iveg_tree_e    = 4
iveg_nontree_s = 5
iveg_nontree_e = 17

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

; _____________________________ Calc Bin Stats ________________________________
          ; 2001, 2017
ts_s    = (/ 366, 6210/)
          ; 2009, 2019
ts_e    = (/3652, 7304/)

; =================================== Path =====================================
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"+case_name(case_num)
end do

; ================================== SET Dim ===================================
filename  = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f         = addfile (filename,"r")

lat_size  = dimsizes(f->y)
lon_size  = dimsizes(f->x)
CABLE_lat = f->y
CABLE_lon = f->x
default_value = f->Rainf@_FillValue
lat_sum = dimsizes(f->Rainf(0,{A:B},0))
lon_sum = dimsizes(f->Rainf(0,0,{C:D}))

Rain     = new((/case_sum, time_tot, lat_sum, lon_sum/), float, default_value)
Evap     = new((/case_sum, time_tot, lat_sum, lon_sum/), float, default_value)
Rain     = 0.
Evap     = 0.

; PFT mask
PFT_mask = where((f->iveg({A:B},{C:D}) .ge. iveg_tree_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_tree_e), 1, f->iveg({A:B},{C:D}))
PFT_mask = where((PFT_mask .ge. iveg_nontree_s) .and. (PFT_mask .le. iveg_nontree_e), 2, PFT_mask)

delete(f)

; ================================= Read Data ==================================
do case_num = 0,case_sum-1
  day_step_s = 0
  do year = year_s,year_e
    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if
    filename1  = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
    f1         = addfile (filename1,"r")
    yr_tot_day = sum(dom)
    day_step_e = day_step_s + yr_tot_day - 1

    Rain(case_num,day_step_s:day_step_e,:,:) = f1->Rainf(:,{A:B},{C:D})
    Evap(case_num,day_step_s:day_step_e,:,:) = f1->Evap(:,{A:B},{C:D})

    day_step_s = day_step_e + 1
    delete(f1)
  end do ; year
end do ; case_num
Rain = Rain*3600*24
Evap = Evap*3600*24

if PFT .eq. "tree" then
  do case_num = 0, case_sum - 1
    do day = 0,time_tot -1
      Rain(case_num,day,:,:)   = where(PFT_mask .eq. 1, Rain(case_num,day,:,:), default_value)
      Evap(case_num,day,:,:)   = where(PFT_mask .eq. 1, Evap(case_num,day,:,:), default_value)
    end do
  end do
else if PFT .eq. "nontree" then
  do case_num = 0, case_sum - 1
    do day = 0,time_tot -1
      Rain(case_num,day,:,:)   = where(PFT_mask .eq. 2, Rain(case_num,day,:,:), default_value)
      Evap(case_num,day,:,:)   = where(PFT_mask .eq. 2, Evap(case_num,day,:,:), default_value)
    end do
  end do
end if
end if

; ================================ READ GLEAM ==================================
; use 2000-2019
if GLEAM_version .eq. "3.3a" then
  path_GLEAM    = "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_3/3_3a/daily/";
  gleam_ET      = path_GLEAM+year_s+"/E_"+year_s+"_GLEAM_v3.3a.nc"
else if GLEAM_version .eq. "3.5a" then
  path_GLEAM    = "/g/data/ua8/GLEAM_v3-5/v3-5a/daily/";
  gleam_ET      = path_GLEAM+year_s+"/E_"+year_s+"_GLEAM_v3.5a.nc"
end if
end if

G_ET          = addfile(gleam_ET,"r")
GLEAM_def_val = doubletofloat(G_ET->E@_FillValue)
tmp           = G_ET->E(:,{C:D},{A:B})
GLEAM_lat     = tmp&lat
GLEAM_lon     = tmp&lon
delete(tmp)
print(GLEAM_def_val)
print(GLEAM_lat)
print(GLEAM_lon)

Evap_GLEAM    = new((/time_tot, lat_sum, lon_sum/), float, GLEAM_def_val)
Evap_GLEAM    = GLEAM_def_val

day_step_s = 0
do year = year_s,year_e
  if ( mod(year,4) .eq. 0) then
     dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
  else
     dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  end if
  yr_tot_day = sum(dom)
  if GLEAM_version .eq. "3.3a" then
    gleam_ET   = path_GLEAM+year+"/E_"+year+"_GLEAM_v3.3a.nc"
  else if GLEAM_version .eq. "3.5a" then
    gleam_ET   = path_GLEAM+year+"/E_"+year+"_GLEAM_v3.5a.nc"
  end if
  end if
  G_ET       = addfile (gleam_ET,"r")
  day_step_e = day_step_s + yr_tot_day - 1
  E_temp     = G_ET->E(:,{C:D}, {A:B})
  E_tmp      = doubletofloat(E_temp(time|:,lat|:, lon|:))
  delete(E_temp)

  ; Regrid
  Evap_GLEAM(day_step_s:day_step_e,:,:) = \
  linint2_Wrap(GLEAM_lon,GLEAM_lat,E_tmp,False,CABLE_lon,CABLE_lat,0)

  delete(gleam_ET)
  delete(G_ET)
  delete(E_tmp)
  day_step_s = day_step_e + 1
end do
Evap_GLEAM = Evap_GLEAM

; ___________________________ Mask GLEAM ______________________________
if PFT .eq. "tree" then
  do day = 0,time_tot -1
    Evap_GLEAM(day,:,:) = where(PFT_mask .eq. 1 .and. .not. ismissing(Evap_GLEAM(day,:,:)),\
                                Evap_GLEAM(day,:,:), default_value)
  end do
else if PFT .eq. "nontree" then
  do day = 0,time_tot -1
    Evap_GLEAM(day,:,:) = where(PFT_mask .eq. 2 .and. .not. ismissing(Evap_GLEAM(day,:,:)),\
                                Evap_GLEAM(day,:,:), default_value)
  end do
end if
end if
Evap_GLEAM@_FillValue = default_value

; ============================= plot time series ===============================
total = 0
do year = 2001,2009
  if mod(year,4) .eq. 0 then
    total = total + 366
  else
    total = total + 365
  end if
end do

P     = dim_avg_n_Wrap(Rain,(/2,3/))
E     = dim_avg_n_Wrap(Evap,(/2,3/))
GLEAM_E = dim_avg_n_Wrap(Evap_GLEAM,(/1,2/))
P_E   = dim_avg_n_Wrap(Rain,(/2,3/))
P_E   = P_E - dim_avg_n_Wrap(Evap,(/2,3/))
P_GLEAM_E = dim_avg_n_Wrap(Rain(0,:,:,:),(/1,2/))
P_GLEAM_E = P_GLEAM_E - dim_avg_n_Wrap(Evap_GLEAM,(/1,2/))
print(P_GLEAM_E)
P_E_acl = new((/(case_sum+1)*2,total/),float, default_value)
P_E_acl = default_value

count = 0
do cnt = 0, dimsizes(ts_s) -1
  do case_num = 0, case_sum -1
    P_E_acl(count,ts_s(cnt)-ts_s(cnt):ts_e(cnt)-ts_s(cnt)) = dim_cumsum_n_Wrap(P_E(case_num,ts_s(cnt):ts_e(cnt)), 0, 0)
    count = count + 1
  end do
  P_E_acl(count,ts_s(cnt)-ts_s(cnt):ts_e(cnt)-ts_s(cnt)) = dim_cumsum_n_Wrap(P_GLEAM_E(ts_s(cnt):ts_e(cnt)), 0, 0)
  count = count + 1
end do

r_CTL_MD = escorc(P(0,ts_s(0):ts_e(0)), E(0,ts_s(0):ts_e(0)))
r_CTL_RD = escorc(P(0,ts_s(1):ts_e(1)), E(0,ts_s(1):ts_e(1)))

r_FD_MD = escorc(P(1,ts_s(0):ts_e(0)), E(1,ts_s(0):ts_e(0)))
r_FD_RD = escorc(P(1,ts_s(1):ts_e(1)), E(1,ts_s(1):ts_e(1)))

r_GLEAM_MD = escorc(P(0,ts_s(0):ts_e(0)), GLEAM_E(ts_s(0):ts_e(0)))
r_GLEAM_RD = escorc(P(0,ts_s(1):ts_e(1)), GLEAM_E(ts_s(1):ts_e(1)))

print("r_CTL_MD: "+r_CTL_MD+"  r_CTL_RD: "+r_CTL_RD)
print("r_FD_MD: "+r_FD_MD+"  r_FD_RD: "+r_FD_RD)
print("r_GLEAM_MD: "+r_GLEAM_MD+"  r_GLEAM_RD: "+r_GLEAM_RD)

; ______________________________ Plotting ______________________________
pic = "./plots/HESS_time_accumulate_P-E_GLEAM_"+GLEAM_version+"_PFT-"+PFT+"_"+year_s+"-"+year_e
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

res                         = True
res@gsnDraw                 = False             ; don't draw yet
res@gsnFrame                = False             ; don't advance frame yet
res@gsnMaximize             = False

res@vpHeightF 	            = 0.4               ; change aspect ratio of plot
res@vpWidthF 	              = 0.8
; res@vpXF                    = 0.1
; res@vpYF                    = 0.9
res@txFontHeightF           = 0.015            ; change title font heights
res@tiYAxisFontThicknessF   = 2.
res@tmXBMajorThicknessF     = 2.
res@tmXBMinorThicknessF     = 2.
res@tmYLMajorThicknessF     = 2.
res@tmYLMinorThicknessF     = 2.
res@tmYLLabelFontThicknessF = 2.
res@tmBorderThicknessF      = 2.
; res@txFont                  = 25

res@lgPerimOn              = False      ; No legend perimeter.
res@pmLegendParallelPosF    = 0.36       ; Move to right
res@pmLegendOrthogonalPosF  = -1.27     ; Move into plot
res@pmLegendWidthF          = 0.24                     ; resize legend width
res@pmLegendHeightF         = 0.21                     ; and height
res@lgBoxMinorExtentF       = 0.28       ; Shorten the legend lines

; res@tiXAxisFontHeightF      = 0.015
; res@tiYAxisFontHeightF      = 0.015
; res@tmXBLabelFontHeightF    = 0.015
; res@tmYLLabelFontHeightF    = 0.015
; res@gsnPaperOrientation     = "portrait";"portrait";"landscape" ;"auto"
res@pmLegendDisplayMode     = "Always"

res@trXMinF                 = 0
res@trXMaxF                 = total
res@trYMinF                 = -110.
res@trYMaxF                 = 800.

res@tmXBMode                = "Explicit"
res@tmXBValues              = ispan(1,total,366)
res@tmXBLabels              = (/"2001~C~2017","","2003~C~2019","","2005","","2007","","2009"/)
x                           = ispan(1,total,1)

res@xyLineThicknesses      = (/ 2., 2., 2.,  2., 2., 2./)
res@xyDashPatterns         = (/  0,  0,  0,   0,  0,  0/)
res@xyLineOpacities        = (/ 1., 1., 1.,  1., 1., 1./) ;0.5, 0.5, 0.5/)
res@xyLineColors           = (/"darkgreen","orange", "black",  "green3","gold","gray50"  /)
res@xyExplicitLegendLabels = (/"GW MD (r="+sprintf("%4.2f",r_CTL_MD)+")","FD MD  (r="+sprintf("%4.2f",r_FD_MD)+")","GLEAM MD", "GW RD (r="+sprintf("%4.2f",r_CTL_RD)+")","FD RD  (r="+sprintf("%4.2f",r_FD_RD)+")", "GLEAM RD"/)  ; create explicit labels
res@lgItemOrder            = (/2,5, 0,3, 1,4/)
res@lgJustification        = "TopRight"


; res@tiXAxisString        = "months"
;**************** ͼ�������� *****************
pres                            = True
pres@wkWidth                    = 2500
pres@wkHeight                   = 2500
pres@gsnMaximize                = False
pres@gsnPanelYWhiteSpacePercent = 0.1
pres@gsnPanelXWhiteSpacePercent = 0.1
; pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
; pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
; pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
; pres@lbTitleFont        = 25
; pres@lbLabelFont        = 25
; pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
; pres@lbTitleFontHeightF = 0.2
; pres@lbLabelFontHeightF = 0.2

plots = new(1,graphic)

res@tiYAxisString = "Accumulated P-E (mm)"
plots             = gsn_csm_xy(wks,x,P_E_acl,res)


; ________________ add order name ________________
txres                  = True
txres@txPerimOn        = False
txres@txFont           = 25
txres@txFontHeightF    = 0.025
amres                  = True
amres@amParallelPosF   = -0.487    ; This is the right edge of the plot.
amres@amOrthogonalPosF = -0.46    ; This is the bottom edge of the plot.
amres@amJust           = "TopLeft"

txid1   = gsn_create_text(wks, "(b)", txres)
annoid1 = gsn_add_annotation(plots(0), txid1, amres)

gsn_panel(wks,(/plots/),(/1,1/),pres)
; delete(plot)
exit
end
