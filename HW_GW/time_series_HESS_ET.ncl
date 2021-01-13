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
time_scale = "daily" ;"monthly" ; "daily"
scale      = "SE-AU"
year_s     = 2000
year_e     = 2018
year_sum   = year_e-year_s+1
time_tot   = year_sum*12 ; plot monthly Evap

; ______________________________ Path _________________________________
path_GLEAM     = "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_3/3_3a/monthly";
path_mum       = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_name      = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr"/)
                  ;"CTL-30x3+30yr-DeepRoot/outputs"/) ; "CTL-30x3+30yr-DeepRoot/outputs-daily"
pic_message    = "CTL_vs_FD"
pic_name       = "HESS_Evap_"+pic_message+"_"+time_scale+"_"+scale+"_"+year_s+"-"+year_e
Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Depth          = (/1.1, 5.1, 15.7, 43.85, 118.55, 316.4, 2355.36/)
                 ; avg(GWdz) over SE Aus = 37.9072 m
                 ; -> last depth=460+(37.9072*10)/2 = 2355.36 (cm)
order_name     = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

; "tree" PFT
iveg_s         = 1
iveg_e         = 4


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
filename1   = path(0)+"/cable_out_1990_SE_Aus.nc"
f1          = addfile (filename1,"r")
PFT_mask    = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
default_val = f1->Qle@_FillValue

delete(filename1)
delete(f1)


; ================================ READ GLEAM ==================================
; use 2000-2019
gleam_ET = path_GLEAM+"/E_1980_2018_GLEAM_v3.3a_MO.nc"
G_ET     = addfile (gleam_ET,"r")
month_s  = (2000-1980)*12

E_tmp    = G_ET->E(month_s:,{C:D},{A:B})

Evap_GLEAM = dim_avg_n_Wrap(E_tmp,(/1,2/))
printVarSummary(Evap_GLEAM)
delete(G_ET)

; ___________________________ READ DATA ______________________________
Evap       = new((/case_sum,2,year_sum*12/),float,default_val)
Evap       = 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  i = 0
  j = 0
  do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
    else
       dom = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
    end if

    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
    f1 = addfile (filename1,"r")

    if time_scale .eq. "monthly" then
      ts = 12
    else if time_scale .eq. "daily" then
      ts = dom(12)
    end if
    end if

    Evap_tmp = f1->Evap(:,{A:B},{C:D})

    do month_num = 0,11
      ; all
      Evap(case_num, 0, j)  = avg(dim_sum_n_Wrap(Evap_tmp(dom(month_num):dom(month_num+1)-1,:,:),0))
      ; trees
      Evap(case_num, 1, j)  = avg(where(PFT_mask .eq. 1, \
                                       dim_sum_n_Wrap(Evap_tmp(dom(month_num):dom(month_num+1)-1,:,:), 0),\
                                       default_val))
      j = j + 1
    end do ; month_num
    delete(Evap_tmp)
    i = i + ts
  end do ; year
end do ; case_num
print("prepare to plot")


Evap_monthly = new((/case_sum+1,year_sum*12/),float,default_val)
Evap_monthly(0,:) = doubletofloat(Evap_GLEAM)
Evap_monthly(1:,:) = Evap(:, 0, :)*24.*3600.


;_____________________________ PLOT ___________________________________
pic = "./plots/"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")

res                         = True
res@gsnDraw                 = False  ; Don't draw plots
res@gsnFrame                = False  ; ����ҳ
res@gsnMaximize             = False

res@vpHeightF               = 0.4
res@vpWidthF                = 0.8
; res@vpXF                    = 0.03
; res@vpYF                    = 0.7
res@tiYAxisFontThicknessF   = 3.
res@tmXBMajorThicknessF     = 3.
res@tmXBMinorThicknessF     = 3.
res@tmYLMajorThicknessF     = 3.
res@tmYLMinorThicknessF     = 3.
res@tmYLLabelFontThicknessF = 3.
res@tmBorderThicknessF      = 3.
res@txFont                  = 25
res@txFontHeightF           = 0.015            ; change title font heights
res@gsnPaperOrientation     = "landscape";"portrait";"landscape" ;"auto"

res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015

res@tmYLMode                = "Explicit"
res@tmXBMode                = "Explicit"
res@txString                = ""
res@gsnLeftString           = ""
res@gsnRightString          = ""

x                           = ispan(1,time_tot,1)
res@tmXBValues              = ispan(1,time_tot,12*2)
res@tmXBLabels              = ispan(year_s,year_e,2)

res@tiYAxisString           = "Evaporation (mm mth~S~-1~N~)"


res@lgPerimOn              = False      ; No legend perimeter.
res@pmLegendParallelPosF   = 0.5       ; Move to right
res@pmLegendOrthogonalPosF = 0.5     ; Move into plot
res@pmLegendWidthF         = 0.3       ; Decrease width
res@pmLegendHeightF        = 0.1       ; Decrease height
res@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines
res@lgItemOrder            = (/ 2,1,0 /)   ; Reorder the legends
res@xyLineColors           = (/"black","blue","red"/)
res@xyExplicitLegendLabels = (/"GLEAM","CTL","FD"/)
plot                       = gsn_csm_xy(wks, x, Evap_monthly,res)

txres                      = True
txres@txPerimOn            = False
txres@txFont               = 25
txres@txFontHeightF        = 0.03

amres                      = True
amres@amParallelPosF       = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF     = -0.46    ; This is the bottom edge of the plot.
amres@amJust               = "TopLeft"

txid1   = gsn_create_text(wks, order_name(0), txres)
annoid1 = gsn_add_annotation(plot, txid1, amres)
draw(wks)
frame(wks)

;**************** ͼ�������� *****************
; pres                    = True                                      ; ���� panel resources.
; pres@gsnMaximize        = True
; pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
; pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
; pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
; pres@lbTitleFont        = 25
; pres@lbLabelFont        = 25
; pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
; pres@lbTitleFontHeightF = 0.2
; pres@lbLabelFontHeightF = 0.2


; plots              = new(6,graphic)

end
