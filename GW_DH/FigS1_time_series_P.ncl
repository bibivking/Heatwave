
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
case_name = "CTL-30x3+30yr/outputs-daily"
time_scale= "monthly"
PFT       = "all" ; "tree" ;"all"
scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
case_sum  = dimsizes(case_name)
order_name  = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

year_s   = 2000
year_e   = 2019
year_sum = year_e-year_s+1
time_tot = year_sum*12

; set PFT
iveg_tree_s    = 1
iveg_tree_e    = 4
iveg_nontree_s = 5
iveg_nontree_e = 17

; =================================== Path =====================================
path = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"+case_name


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
ts_e    = (2019-year_s)*12+11 ;(/ (2009-year_s)*12+11, (2019-year_s)*12+11/)

; ================================== SET Dim ===================================
filename  = path + "/cable_out_"+year_s+"_SE_Aus.nc"
f         = addfile (filename,"r")

lat_size  = dimsizes(f->y)
lon_size  = dimsizes(f->x)
CABLE_lat = f->y
CABLE_lon = f->x
default_value = f->Rainf@_FillValue
lat_sum = dimsizes(f->Rainf(0,{A:B},0))
lon_sum = dimsizes(f->Rainf(0,0,{C:D}))

Rain     = new((/time_tot,lat_sum,lon_sum/), float, default_value)
Rain     = 0.

; PFT mask
PFT_mask = where((f->iveg({A:B},{C:D}) .ge. iveg_tree_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_tree_e), 1, f->iveg({A:B},{C:D}))
PFT_mask = where((PFT_mask .ge. iveg_nontree_s) .and. (PFT_mask .le. iveg_nontree_e), 2, PFT_mask)

delete(f)


; ================================= Read Data ==================================
cnt = 0
do year = year_s,year_e
  print("year = " +year)
  if ( mod(year,4) .eq. 0) then
     dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
  else
     dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  end if
  filename1  = path+"/cable_out_"+year+"_SE_Aus.nc"
  f1         = addfile (filename1,"r")

  if year .lt. 2000 then
    tmp = f1->Rainf(:,{A:B},{C:D})
    do month = 0,11
      Rain(cnt,:,:) = tmp(month,:,:)*dom(month)
      cnt = cnt + 1
    end do
    delete(tmp)
  else
    tmp = f1->Rainf(:,{A:B},{C:D})
    day_step_s = 0
    do month = 0,11
      day_step_e = sum(dom(0:month))-1
      Rain(cnt,:,:) = dim_sum_n_Wrap(tmp(day_step_s:day_step_e,:,:),0)
      day_step_s = day_step_e + 1
      cnt = cnt + 1
    end do
    delete(tmp)
  end if
  delete(f1)
end do ; year
Rain = Rain*3600*24

Rain_avg = dim_avg_n_Wrap(Rain,(/1,2/))

; ================================== Check =====================================
pic1 = "./plots/FigS1_Check"+year_s+"-"+year_e
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

res1                = True
res1@gsnDraw        = True             ; don't draw yet
res1@gsnFrame       = True             ; don't advance frame yet
res1@gsnMaximize    = False

res1@cnFillMode     = "RasterFill"            ; Raster Mode
res1@cnFillOn       = True                            ;��ɫ
res1@cnLinesOn      = False
res1@cnLineLabelsOn = False

plots = gsn_csm_contour(wks1, Rain(0,:,:), res1)


; ================================== Plotting ==================================
pic = "./plots/FigS1_time_series_P_"+year_s+"-"+year_e
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

res                         = True
res@gsnDraw                 = False             ; don't draw yet
res@gsnFrame                = False             ; don't advance frame yet
res@gsnMaximize             = True

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
res@pmLegendOrthogonalPosF  = -1.22 ;7     ; Move into plot
res@pmLegendWidthF          = 0.1                     ; resize legend width
res@pmLegendHeightF         = 0.1                     ; and height
res@lgBoxMinorExtentF       = 0.3       ; Shorten the legend lines
res@pmLegendDisplayMode     = "Never"

res@trXMinF                 = ts_s
res@trXMaxF                 = ts_e+1

res@tmXBMode                = "Explicit"
res@tmXBValues              = ispan(ts_s,ts_e+1,24);ispan(ts_s,ts_e,24)
res@tmXBLabels              = (/"2000","2002","2004","2006","2008","2010","2012","2014","2016","2018","2020"/)
x                           = ispan(ts_s,ts_e,1);ispan(ts_s,ts_e,1)

res@xyLineThicknessF      = 1.5
res@xyDashPattern         = 0
res@xyLineOpacityF        = 1.
res@xyLineColor           = "black"
res@gsnXYBarChartColors   = "gray80"
; res@xyExplicitLegendLabels = (/"All","Shallow","Medium","Deep"/)
; res@lgJustification        = "TopRight"


;**************** ͼ�������� *****************
pres                            = True
pres@wkWidth                    = 2500
pres@wkHeight                   = 2500
pres@gsnMaximize                = False
pres@gsnPanelYWhiteSpacePercent = 0.1
pres@gsnPanelXWhiteSpacePercent = 0.1

plots = new(1,graphic)

res@trYMinF                = 0
res@trYMaxF                = 130
res@gsnXYBarChart          = True            ; create bar chart
res@tiYAxisString          = "P (mm mth~S~-1~N~)"
plots(0)                   = gsn_csm_xy(wks,x,Rain_avg,res)
print(x)
; ________________ add shading ________________
poly             = new(2,graphic)
gres             = True;"rosybrown" ;"azure"
gres@gsEdgesOn   = False
gres@gsEdgeColor = "red"
gres@gsEdgeDashPattern = 0
gres@gsEdgeThicknessF  = 2
gres@gsFillColor = "gray80";"cornsilk"

xmin = (/(2001-year_s)*12,(2017-year_s)*12/)
xmax = (/(2010-year_s)*12,(2020-year_s)*12/)
ymin = (/  0,   0/)
ymax = (/130, 130/)

xbars = (/ (/xmin(0), xmax(0), xmax(0), xmin(0), xmin(0)/), \
            (/xmin(1), xmax(1), xmax(1), xmin(1), xmin(1)/) /)
print(xbars)
ybars = (/ (/ymin(0),ymin(0),ymax(0),ymax(0),ymin(0)/), \
           (/ymin(1),ymin(1),ymax(1),ymax(1),ymin(1)/)  /)
print(ybars)
counter = 0
do shadow_num = 0,1
  poly(counter) = gsn_add_polygon(wks,plots(0),xbars(shadow_num,:),ybars(shadow_num,:),gres)
  counter = counter + 1
end do

gsn_panel(wks,(/plots/),(/1,1/),pres)
exit
end
