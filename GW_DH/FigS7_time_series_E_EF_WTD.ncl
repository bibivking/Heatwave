
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
case_name = (/"CTL-30x3+30yr/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr"/)

scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
case_sum  = dimsizes(case_name)
order_name  = (/"(a)","(b)","(c)","(d)","(e)","(f)"/) ;

year_s   = 2000
year_e   = 2019
year_sum = year_e-year_s+1
time_tot = year_sum*12

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
lat_sum = dimsizes(f->Rainf(0,{A:B},0))
lon_sum = dimsizes(f->Rainf(0,0,{C:D}))
default_value = f->Rainf@_FillValue
delete(f)

; ================================= Read Data ==================================
Evap     = new((/case_sum, time_tot, lat_sum, lon_sum/), float, default_value)
Qle      = new((/case_sum, time_tot, lat_sum, lon_sum/), float, default_value)
Qh       = new((/case_sum, time_tot, lat_sum, lon_sum/), float, default_value)
WTD      = new((/case_sum, time_tot, lat_sum, lon_sum/), float, default_value)
Evap     = 0.
Qle      = 0.
Qh       = 0.
WTD      = 0.

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
      Evap_tmp = f1->Evap(:,{A:B},{C:D})
      Qle_tmp  = f1->Qle(:,{A:B},{C:D})
      Qh_tmp   = f1->Qh(:,{A:B},{C:D})
      WTD_tmp  = f1->WatTable(:,{A:B},{C:D})
      do month = 0,11
        Evap(case_num,cnt,:,:) = Evap_tmp(month,:,:)*dom(month)
        Qle(case_num,cnt,:,:)  = Qle_tmp(month,:,:)
        Qh(case_num,cnt,:,:)   = Qh_tmp(month,:,:)
        WTD(case_num,cnt,:,:)  = WTD_tmp(month,:,:)
        cnt = cnt + 1
      end do
      delete(Evap_tmp)
      delete(Qle_tmp)
      delete(Qh_tmp)
      delete(WTD_tmp)
    else
      Evap_tmp = f1->Evap(:,{A:B},{C:D})
      Qle_tmp  = f1->Qle(:,{A:B},{C:D})
      Qh_tmp   = f1->Qh(:,{A:B},{C:D})
      WTD_tmp  = f1->WatTable(:,{A:B},{C:D})
      day_step_s = 0
      do month = 0,11
        day_step_e = sum(dom(0:month))-1
        Evap(case_num,cnt,:,:) = dim_sum_n_Wrap(Evap_tmp(day_step_s:day_step_e,:,:),0)
        Qle(case_num,cnt,:,:)  = dim_avg_n_Wrap(Qle_tmp(day_step_s:day_step_e,:,:),0)
        Qh(case_num,cnt,:,:)   = dim_avg_n_Wrap(Qh_tmp(day_step_s:day_step_e,:,:),0)
        WTD(case_num,cnt,:,:)  = dim_avg_n_Wrap(WTD_tmp(day_step_s:day_step_e,:,:),0)
        day_step_s = day_step_e + 1
        cnt = cnt + 1
      end do
      delete(Evap_tmp)
      delete(Qle_tmp)
      delete(Qh_tmp)
      delete(WTD_tmp)
    end if
    delete(f1)
  end do ; year
end do ; case_num
Evap = Evap*3600*24

; ================================== Check =====================================
pic1 = "./plots/FigS7_Check_"+year_s+"-"+year_e
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

res1                = True
res1@gsnDraw        = True             ; don't draw yet
res1@gsnFrame       = True             ; don't advance frame yet
res1@gsnMaximize    = False

res1@cnFillOn                = True
res1@cnMissingValFillColor   = "gray80"
res1@cnFillBackgroundColor   = "gray80"

res1@cnFillMode     = "RasterFill"            ; Raster Mode
res1@cnFillOn       = True                            ;��ɫ
res1@cnLinesOn      = False
res1@cnLineLabelsOn = False
res1@cnLevelSelectionMode = "ExplicitLevels"

res1@cnLevels        = fspan(0,200,11)
plot = gsn_csm_contour(wks1, Evap(0,6,:,:), res1)
delete(res1@cnLevels)

res1@cnLevels        = fspan(0,1,11)
plot = gsn_csm_contour(wks1, Qle(0,6,:,:), res1)
delete(res1@cnLevels)

res1@cnLevels        = fspan(0,1,11)
plot = gsn_csm_contour(wks1, Qh(0,6,:,:), res1)
delete(res1@cnLevels)

res1@cnLevels        = fspan(0,30,16)
plot = gsn_csm_contour(wks1, WTD(0,6,:,:), res1)
delete(res1@cnLevels)

; ================================== Mask ======================================
WTD_mask       = dim_avg_n_Wrap(WTD(0,:,:,:),0) ; GW WTD
WTD_mask_multi = new((/time_tot,lat_sum,lon_sum/), float, default_value)
do z = 0, time_tot-1
  WTD_mask_multi(z,:,:) = WTD_mask
end do
delete(WTD_mask)

WTD_lvl      = 4 ;"all", "shallow"; "medium"; "deep"

deltaE_class = new((/WTD_lvl, time_tot, lat_sum, lon_sum/), float, default_value)
GW_Qle_class = new((/WTD_lvl, time_tot, lat_sum, lon_sum/), float, default_value)
GW_Qh_class  = new((/WTD_lvl, time_tot, lat_sum, lon_sum/), float, default_value)
GW_WTD_class = new((/WTD_lvl, time_tot, lat_sum, lon_sum/), float, default_value)

; deltaE
deltaE_class(0,:,:,:) = Evap(0,:,:,:) - Evap(1,:,:,:)
deltaE_class(1,:,:,:) = where( WTD_mask_multi .le. 5., Evap(0,:,:,:) - Evap(1,:,:,:), default_value )
deltaE_class(2,:,:,:) = where( WTD_mask_multi .gt. 5. .and. WTD_mask_multi .le. 10.,\
                               Evap(0,:,:,:) - Evap(1,:,:,:), default_value )
deltaE_class(3,:,:,:) = where( WTD_mask_multi .gt. 10., Evap(0,:,:,:) - Evap(1,:,:,:), default_value )

; GW Qle
GW_Qle_class(0,:,:,:) = Qle(0,:,:,:)
GW_Qle_class(1,:,:,:) = where( WTD_mask_multi .le. 5., Qle(0,:,:,:), default_value )
GW_Qle_class(2,:,:,:) = where( WTD_mask_multi .gt. 5. .and. WTD_mask_multi .le. 10.,\
                              Qle(0,:,:,:), default_value )
GW_Qle_class(3,:,:,:) = where( WTD_mask_multi .gt. 10., Qle(0,:,:,:), default_value )

; GW Qh
GW_Qh_class(0,:,:,:) = Qh(0,:,:,:)
GW_Qh_class(1,:,:,:) = where( WTD_mask_multi .le. 5., Qh(0,:,:,:), default_value )
GW_Qh_class(2,:,:,:) = where( WTD_mask_multi .gt. 5. .and. WTD_mask_multi .le. 10.,\
                              Qh(0,:,:,:), default_value )
GW_Qh_class(3,:,:,:) = where( WTD_mask_multi .gt. 10., Qh(0,:,:,:), default_value )

; GW WTD
GW_WTD_class(0,:,:,:) = WTD(0,:,:,:)
GW_WTD_class(1,:,:,:) = where( WTD_mask_multi .le. 5., WTD(0,:,:,:), default_value )
GW_WTD_class(2,:,:,:) = where( WTD_mask_multi .gt. 5. .and. WTD_mask_multi .le. 10.,\
                              WTD(0,:,:,:), default_value )
GW_WTD_class(3,:,:,:) = where( WTD_mask_multi .gt. 10., WTD(0,:,:,:), default_value )

delete(WTD_mask_multi)

; Check
res1@cnLevels        = fspan(0,200,11)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(deltaE_class(1,:11,:,:),0), res1)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(deltaE_class(2,:11,:,:),0), res1)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(deltaE_class(3,:11,:,:),0), res1)
delete(res1@cnLevels)
; res1@cnLevels        = fspan(0,1,11)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(GW_Qle_class(1,:11,:,:),0), res1)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(GW_Qle_class(2,:11,:,:),0), res1)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(GW_Qle_class(3,:11,:,:),0), res1)
; res1@cnLevels        = fspan(0,1,11)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(GW_Qh_class(1,:11,:,:),0), res1)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(GW_Qh_class(2,:11,:,:),0), res1)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(GW_Qh_class(3,:11,:,:),0), res1)
delete(res1@cnLevels)
res1@cnLevels        = fspan(0,30,16)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(GW_WTD_class(1,:11,:,:),0), res1)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(GW_WTD_class(2,:11,:,:),0), res1)
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(GW_WTD_class(3,:11,:,:),0), res1)
delete(res1@cnLevels)

; ======================= regional avg & month2year ============================
; Regional avg
deltaE        = dim_avg_n_Wrap(deltaE_class,(/2,3/))
GW_Qle        = dim_avg_n_Wrap(GW_Qle_class,(/2,3/))
GW_Qh         = dim_avg_n_Wrap(GW_Qh_class,(/2,3/))
GW_WTD        = dim_avg_n_Wrap(GW_WTD_class,(/2,3/))
print(GW_Qle)
print(GW_Qh)

; monthly to annual
deltaE_yr     = new((/WTD_lvl, year_sum/), float, default_value)
GW_EF_yr      = new((/WTD_lvl, year_sum/), float, default_value)
GW_WTD_yr     = new((/WTD_lvl, year_sum/), float, default_value)

do lvl = 0,3
  deltaE_yr(lvl,:) = month_to_annual(deltaE(lvl,:),0)
  GW_EF_yr(lvl,:)  = month_to_annual(GW_Qle(lvl,:)/(GW_Qle(lvl,:)+GW_Qh(lvl,:)),1)
  GW_WTD_yr(lvl,:) = month_to_annual(GW_WTD(lvl,:),1)
end do
print(GW_EF_yr)
; ______________________________ Plotting ______________________________
pic = "./plots/FigS7_time_series_E_EF_WTD_"+year_s+"-"+year_e
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

res                         = True
res@gsnDraw                 = True             ; don't draw yet
res@gsnFrame                = True             ; don't advance frame yet
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


; res@tiXAxisFontHeightF      = 0.015
; res@tiYAxisFontHeightF      = 0.015
; res@tmXBLabelFontHeightF    = 0.015
; res@tmYLLabelFontHeightF    = 0.015
; res@gsnPaperOrientation     = "portrait";"portrait";"landscape" ;"auto"

res@trXMinF                 = year_s
res@trXMaxF                 = year_e

res@tmXBMode                = "Explicit"
res@tmXBValues              = ispan(year_s,year_e,2)
res@tmXBLabels              = (/"2000","2002","2004","2006","2008","2010","2012","2014","2016","2018"/)
x                           = ispan(year_s,year_e,1)

res@xyLineThicknesses      = (/  2,  2,  2,  2   /)
res@xyDashPatterns         = (/  0,  0,  0,  0   /)
res@xyLineOpacities        = (/  1,  1,  1,  1   /)
res@xyLineColors           = (/"black","blue","darkgreen","orange"/)
res@xyExplicitLegendLabels = (/"All","Shallow","Medium","Deep"/)
res@lgItemOrder            = (/3,2,1,0/)
res@lgJustification        = "TopRight"

res@lgPerimOn              = False      ; No legend perimeter.
res@pmLegendParallelPosF    = 0.55       ; Move to right
res@pmLegendOrthogonalPosF  =-1.25         ; Move into plot
res@pmLegendWidthF          = 0.22                     ; resize legend width
res@pmLegendHeightF         = 0.14                     ; and height
res@lgBoxMinorExtentF       = 0.3       ; Shorten the legend lines


; res@tiXAxisString        = "months"
;**************** ͼ�������� *****************
pres                            = True
pres@wkWidth                    = 2500
pres@wkHeight                   = 2500
pres@gsnMaximize                = False
pres@gsnPanelYWhiteSpacePercent = 0.1
pres@gsnPanelXWhiteSpacePercent = 0.1

plots = new(3,graphic)

res@pmLegendDisplayMode  = "Always"
res@gsnXYBarChart        = False

res@tiYAxisString = "~F5~D"+"~F25~E (mm yr~S~-1~N~)"
res@trYMinF       = -10.
res@trYMaxF       = 150.
plots(0)          = gsn_csm_xy(wks,x,deltaE_yr,res)

res@pmLegendDisplayMode  = "Never"

res@tiYAxisString = "~F25~EF (-)"
res@trYMinF       = 0.
res@trYMaxF       = 0.65
plots(1)          = gsn_csm_xy(wks,x,GW_EF_yr,res)

res@trYReverse    = True
res@tiYAxisString = "~F25~WTD (m)"
res@trYMinF       = 0.
res@trYMaxF       = 13.5
plots(2)          = gsn_csm_xy(wks,x,GW_WTD_yr,res)

; ________________ add order name ________________
txres                  = True
txres@txPerimOn        = False
txres@txFont           = 25
txres@txFontHeightF    = 0.03
amres                  = True
amres@amParallelPosF   = -0.487    ; This is the right edge of the plot.
amres@amOrthogonalPosF = -0.46    ; This is the bottom edge of the plot.
amres@amJust           = "TopLeft"

; do cnt = 0,2
;   txid1   = gsn_create_text(wks, order_name(cnt), txres)
;   annoid1 = gsn_add_annotation(plots(cnt), txid1, amres)
; end do

gsn_panel(wks,(/plots(0)/),(/1,1/),pres)

exit
end
