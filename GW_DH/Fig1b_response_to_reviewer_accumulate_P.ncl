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

year_s   = 2001
year_e   = 2009
year_sum = year_e-year_s+1
time_tot = year_sum*12

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

if time_scale .eq. "monthly" then
            ; 2001, 2017
  ts_s    = (2000-year_s)*12    ;(/ (2001-year_s)*12, (2017-year_s)*12/)
            ; 2009, 2019
  ts_e    = (2019-year_s)*12+11 ;(/ (2009-year_s)*12+11, (2019-year_s)*12+11/)
end if


; ================================== SET Dim ===================================
filename  = path + "/cable_out_"+year_s+"_SE_Aus.nc"
f         = addfile (filename,"r")

lat_size  = dimsizes(f->y)
lon_size  = dimsizes(f->x)
CABLE_lat = f->y
CABLE_lon = f->x
default_value = f->Rainf@_FillValue
print(default_value)
lat_sum = dimsizes(f->Rainf(0,{A:B},0))
lon_sum = dimsizes(f->Rainf(0,0,{C:D}))

Rain     = new((/time_tot,lat_sum,lon_sum/), float, default_value)
Rain     = default_value

delete(f)

; ================================== Check =====================================
pic1 = "./plots/Fig1b_response_to_reviewer_Check_"+year_s+"-"+year_e
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

res1                = True

res1@cnFillOn                = True
res1@cnMissingValFillColor   = "gray80"
res1@cnFillBackgroundColor   = "gray80"
res1@gsnDraw        = True             ; don't draw yet
res1@gsnFrame       = True             ; don't advance frame yet
res1@gsnMaximize    = False

res1@cnFillMode     = "RasterFill"            ; Raster Mode
res1@cnFillOn       = True                            ;��ɫ
res1@cnLinesOn      = False
res1@cnLineLabelsOn = False
plot = gsn_csm_contour(wks1, dim_avg_n_Wrap(Rain(:,:,:),0), res1)

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

  tmp = f1->Rainf
  day_step_s = 0
  do month = 0,11
    day_step_e = sum(dom(0:month))-1
    Rain(cnt,:,:) = dim_sum_n_Wrap(tmp(day_step_s:day_step_e,:,:),0)
    day_step_s = day_step_e + 1
    cnt = cnt + 1
  end do
  delete(tmp)
end do ; year
Rain = Rain*3600*24
Rain@_FillValue = default_value
plot = gsn_csm_contour(wks1, dim_sum_n_Wrap(Rain(:,:,:),0), res1)
Rain_avg = dim_avg_n_Wrap(Rain,(/1,2/))
Rain_avg@_FillValue = default_value
P_acl = dim_cumsum_n_Wrap(Rain_avg, 0, 0)

; __________________________ WRITE TXT OUTPUT _____________________________
option        = True
option@row    = False

option@fout =  "./txt/Fig1b_response_to_reviewer_accumulate_P_"+year_s+"-"+year_e+".txt"
format   = dimsizes(P_acl)+"f15.2"
print(P_acl)
write_matrix(P_acl, format, option)
delete(option@fout)

; ================================== Plotting ==================================
pic = "./plots/Fig1b_response_to_reviewer_accumulate_P_"+time_scale+"_"+year_s+"-"+year_e
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
res@tmXBLabels              = (/"2001","2002","2003","2004","2005","2006","2007","2008","2009"/)
x                           = ispan(ts_s,ts_e,1);ispan(ts_s,ts_e,1)
print(x)



;**************** ͼ�������� *****************
pres                            = True
pres@wkWidth                    = 2500
pres@wkHeight                   = 2500
pres@gsnMaximize                = False
pres@gsnPanelYWhiteSpacePercent = 0.1
pres@gsnPanelXWhiteSpacePercent = 0.1

plots = new(1,graphic)
res@gsnXYBarChart          = True            ; create bar chart

if time_scale .eq. "monthly" then
  res@xyLineThicknessF      = 1.5
  res@xyDashPattern         = 0
  res@xyLineOpacityF        = 1.
  res@xyLineColor           = "black"
  res@gsnXYBarChartColors   = "gray80"
  res@tmXBValues            = ispan(ts_s,ts_e+1,12);ispan(ts_s,ts_e,24)
  res@trYMinF               = 0
  res@trYMaxF               = 130
  res@tiYAxisString         = "accumulate P (mm)"
  plots(0)                  = gsn_csm_xy(wks,x,P_acl,res)
end if

gsn_panel(wks,(/plots/),(/1,1/),pres)
exit
end
