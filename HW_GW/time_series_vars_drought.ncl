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
pic_name  = "droughts"
case_name = (/"CTL-30x3+30yr/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr"/)

PFT       = "all" ; "tree" ;"all"
scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
case_sum  = dimsizes(case_name)
order_name  = (/"(a)","(b)","(c)","(d)","(e)","(f)"/) ; "(a)","(b)",
var_name  = (/"Rainf","Evap","GPP","Qle","Qh","WatTable"/)
WTD       = "shallow" ; "medium", "deep"

year_s   = 2000
year_e   = 2000;19
year_sum = year_e-year_s+1
time_tot = year_sum*12
var_sum  = dimsizes(var_name)

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
ts_s    = (2000-year_s)*12    ;(/ (2001-year_s)*12, (2017-year_s)*12/)
          ; 2009, 2019
ts_e    = (2000-year_s)*12+11 ;(/ (2009-year_s)*12+11, (2019-year_s)*12+11/)
;2019
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

Vars     = new((/case_sum, var_sum, time_tot,lat_sum,lon_sum/), float, default_value)
Vars     = 0.

; PFT mask
PFT_mask = where((f->iveg({A:B},{C:D}) .ge. iveg_tree_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_tree_e), 1, f->iveg({A:B},{C:D}))
PFT_mask = where((PFT_mask .ge. iveg_nontree_s) .and. (PFT_mask .le. iveg_nontree_e), 2, PFT_mask)

delete(f)

; ================================= Read Data ==================================
do case_num = 0,case_sum-1
  cnt = 0
  do year = year_s,year_e
    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if
    filename1  = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
    f1         = addfile (filename1,"r")

    if year .lt. 2000 then
      tmp = new((/var_sum, 12, lat_sum, lon_sum/), float, default_value)
      do var_num = 0, var_sum-1
        tmp(var_num,:,:,:) = f1->$var_name(var_num)$(:,{A:B},{C:D})
      end do
      do month = 0,11
        Vars(case_num,0:2,cnt,:,:) = tmp(0:2,month,:,:)*dom(month)
        Vars(case_num,3:5,cnt,:,:) = tmp(3:5,month,:,:)
        cnt = cnt + 1
      end do
      delete(tmp)
    else
      tmp = new((/var_sum, sum(dom), lat_sum, lon_sum/), float, default_value)
      do var_num = 0, var_sum-1
        tmp(var_num,:,:,:) = f1->$var_name(var_num)$(:,{A:B},{C:D})
      end do
      day_step_s = 0
      do month = 0,11
        day_step_e = sum(dom(0:month))-1
        Vars(case_num,0:2,cnt,:,:) = dim_sum_n_Wrap(tmp(0:2,day_step_s:day_step_e,:,:),1)
        Vars(case_num,3:5,cnt,:,:) = dim_avg_n_Wrap(tmp(3:5,day_step_s:day_step_e,:,:),1)
        day_step_s = day_step_e + 1
        cnt = cnt + 1
      end do
      delete(tmp)
    end if
    delete(f1)
  end do ; year
end do ; case_num
Vars(:,0:2,:,:,:) = Vars(:,0:2,:,:,:)*3600*24

; ================================== Mask ======================================
WTD_mask = dim_avg_n_Wrap(Vars(0,5,:,:,:),0)

WTD_mask_multi = new((/time_tot,lat_sum,lon_sum/), float, default_value)
PFT_mask_multi = new((/time_tot,lat_sum,lon_sum/), float, default_value)
do z = 0, time_tot-1
  WTD_mask_multi(z,:,:) = WTD_mask
  PFT_mask_multi(z,:,:) = PFT_mask
end do
delete(PFT_mask)
delete(WTD_mask)

do case_num = 0, case_sum - 1
  do var_num = 0,var_sum -1
    if PFT .eq. "tree" then
      Vars(case_num,var_num,:,:,:)   = where(PFT_mask_multi .eq. 1, Vars(case_num,var_num,:,:,:), default_value)
    end if
    if PFT .eq. "nontree" then
      Vars(case_num,var_num,:,:,:)   = where(PFT_mask_multi .eq. 2, Vars(case_num,var_num,:,:,:), default_value)
    end if
  end do
end do
delete(PFT_mask_multi)

; calculate diff
Vars_diff = new((/4,var_num,time_tot,lat_sum,lon_sum/), float, default_value)
; Rain
Vars_diff(0,0,:,:,:)   = Vars(0,0,:,:,:)
; E GPP Qle Qh
Vars_diff(0,1:4,:,:,:) = Vars(0,1:4,:,:,:)  - Vars(1,1:4,:,:,:)
; WTD
Vars_diff(0,5,:,:,:)   = Vars(0,5,:,:,:)

do case_num = 0, case_sum - 1
  do var_num = 0,var_sum -1
    ; WTD .eq. "shallow"
    Vars_diff(1,var_num,:,:,:)  = where( WTD_mask_multi .le. 5., Vars(case_num,var_num,:,:,:), default_value )
    ; WTD .eq. "medium"
    Vars_diff(2,var_num,:,:,:)  = where( WTD_mask_multi .gt. 5. .and. WTD_mask_multi .le. 10., Vars(case_num,var_num,:,:,:), default_value )
    ; WTD .eq. "deep"
    Vars_diff(3,var_num,:,:,:)  = where( WTD_mask_multi .gt. 10., Vars(case_num,var_num,:,:,:), default_value )
  end do
end do
delete(WTD_mask_multi)

; ============================= plot time series ===============================
V_diff   = dim_avg_n_Wrap(Vars_diff,(/3,4/))

; deseasonal
; V_deseason   = V ; rainfall doesn't need decycle
; V_month_mean = clmMonLLT(V)

; do year = year_s, year_e
;   month_s = (year - year_s)*12
;   V_deseason(0,1:4,month_s:month_s+11) = V_deseason(0,1:4,month_s:month_s+11) - V_month_mean(0,1:4,:)
;   V_deseason(1,1:4,month_s:month_s+11) = V_deseason(1,1:4,month_s:month_s+11) - V_month_mean(0,1:4,:)
; end do

; ______________________________ Plotting ______________________________
pic = "./plots/HESS_time_series_droughts_PFT-"+PFT+"_"+year_s+"-"+year_e
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
res@pmLegendParallelPosF    = 0.38       ; Move to right
res@pmLegendOrthogonalPosF  = -1.27     ; Move into plot
res@pmLegendWidthF          = 0.1                     ; resize legend width
res@pmLegendHeightF         = 0.1                     ; and height
res@lgBoxMinorExtentF       = 0.3       ; Shorten the legend lines

; res@tiXAxisFontHeightF      = 0.015
; res@tiYAxisFontHeightF      = 0.015
; res@tmXBLabelFontHeightF    = 0.015
; res@tmYLLabelFontHeightF    = 0.015
; res@gsnPaperOrientation     = "portrait";"portrait";"landscape" ;"auto"
res@pmLegendDisplayMode     = "Always"

res@trXMinF                 = ts_s
res@trXMaxF                 = ts_e
; res@trYMinF                 = -110.
; res@trYMaxF                 = 800.

res@tmXBMode                = "Explicit"
res@tmXBValues              = ispan(ts_s,ts_e,24)
res@tmXBLabels              = (/"2000","2002","2004","2006","2008","2010","2012","2014","2016","2018"/)
x                           = ispan(ts_s,ts_e,1)

res@xyLineThicknesses      = (/ 1.5, 1.5, 1.5, 1.5/)
res@xyDashPatterns         = (/  0,  0,  0,  0   /)
res@xyLineOpacities        = (/0.5, 0.5, 0.5, 0.5 /)
res@xyLineColors           = (/"black","darkgreen","orange","red"/)
res@xyExplicitLegendLabels = (/"All","Shallow","Medium","Deep"/)
;(/"GW","FD"/)
; res@lgItemOrder            = (/2,5, 0,3, 1,4/)
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

plots = new(6,graphic)

res@gsnXYBarChart = False            ; create bar chart
res@tiYAxisString = "P (mm/mth)"
plots(0)          = gsn_csm_xy(wks,x,V_diff(:,0,ts_s:ts_e),res)

res@gsnXYBarChart = False
res@tiYAxisString = "E_diff (mm/mth)"
plots(1)          = gsn_csm_xy(wks,x,V_diff(:,1,ts_s:ts_e),res)

res@tiYAxisString = "GPP_diff (umol/m2)"
plots(2)          = gsn_csm_xy(wks,x,V_diff(:,2,ts_s:ts_e),res)

res@tiYAxisString = "Qle_diff (W/m2)"
plots(3)          = gsn_csm_xy(wks,x,V_diff(:,3,ts_s:ts_e),res)

res@tiYAxisString = "Qh_diff (W/m2)"
plots(4)          = gsn_csm_xy(wks,x,V_diff(:,4,ts_s:ts_e),res)

res@tiYAxisString = "WTD_diff (m)"
plots(5)          = gsn_csm_xy(wks,x,V_diff(:,5,ts_s:ts_e),res)

; ________________ add order name ________________
txres                  = True
txres@txPerimOn        = False
txres@txFont           = 25
txres@txFontHeightF    = 0.025
amres                  = True
amres@amParallelPosF   = -0.487    ; This is the right edge of the plot.
amres@amOrthogonalPosF = -0.46    ; This is the bottom edge of the plot.
amres@amJust           = "TopLeft"
do var_num = 0,var_sum-1
  txid1   = gsn_create_text(wks, order_name(var_num), txres)
  annoid1 = gsn_add_annotation(plots(var_num), txid1, amres)
end do
gsn_panel(wks,(/plots/),(/3,2/),pres)
; delete(plot)
exit
end

; ; ================================ READ DOLCE ==================================
; LH_to_ET     = 0.408*0.000001 ; 1 W/m2 = 0.408*0.000001 mm/s, from
;                               ; https://blog.csdn.net/liyanzhong/article/details/46362723
; LH_to_ET_day =  LH_to_ET*3600*24 ; 1 W/m2 = 0.408*0.000001*3600*24 mm/d
;
; path_DOLCE     = "/g/data/w35/Shared_data/Observations/Global_ET_products/DOLCE/v2.1"
; DOLCE_ET       = path_DOLCE+"/DOLCE_v2.1_"+year_s+".nc"
; D_ET           = addfile(DOLCE_ET,"r")
; D_ET_tmp       = D_ET->hfls(0,{A:B},{C:D})
; lat_DOLCE      = D_ET_tmp&lat
; lon_DOLCE      = D_ET_tmp&lon
; Evap_DOLCE_tmp = new((/year_sum*12,dimsizes(lat_DOLCE),dimsizes(lon_DOLCE)/),float,default_value) ;
;
; delete(DOLCE_ET)
; delete(D_ET)
;
; cnt = 0
; do year = year_s,year_e-1
;
;   if ( mod(year,4) .eq. 0) then
;      dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
;   else
;      dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
;   end if
;
;   DOLCE_ET = path_DOLCE+"/DOLCE_v2.1_"+year+".nc"
;   D_ET     = addfile(DOLCE_ET,"r")
;
;   ET_tmp   = D_ET->hfls(:,{A:B},{C:D})
;
;   do month = 0,11
;     Evap_DOLCE_tmp(cnt,:,:) = ET_tmp(month,:,:)*LH_to_ET_day*dom(month)
;     cnt = cnt + 1
;   end do
;
;   delete(ET_tmp)
;   delete(DOLCE_ET)
;   delete(D_ET)
; end do
; replace_ieeenan(Evap_DOLCE_tmp,default_value,0)
;
; ; default value in DOLCE is NAN it may cause issues
; ; ___________________________ Regrid DOLCE ______________________________
; Evap_DOLCE = linint2_Wrap(lon_DOLCE,lat_DOLCE,Evap_DOLCE_tmp,False,CABLE_lon,CABLE_lat,0)
; printVarSummary(Evap_DOLCE)
;
; ; ___________________________ Mask DOLCE ______________________________
; if PFT .eq. "tree" then
;   do month = 0,time_tot -1
;     Evap_DOLCE(month,:,:) = where(PFT_mask .eq. 1 .and. .not. ismissing(Evap_DOLCE(month,:,:)),\
;                                 Evap_DOLCE(month,:,:), Evap_DOLCE@_FillValue)
;   end do
; else if PFT .eq. "nontree" then
;   do month = 0,time_tot -1
;     Evap_DOLCE(month,:,:) = where(PFT_mask .eq. 2 .and. .not. ismissing(Evap_DOLCE(month,:,:)),\
;                                 Evap_DOLCE(month,:,:), Evap_DOLCE@_FillValue)
;   end do
; end if
; end if
; ; Evap_DOLCE@_FillValue = default_value
