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
PFT        = "tree" ; "tree"
scale      = "SE-AU"
year_s     = 2019
year_e     = 2019
year_sum   = year_e-year_s+1
time_tot   = year_sum*365*8
plot_ts_s  = (11-1)*8
plot_ts_e  = (25-1)*8+7
plot_ts    = plot_ts_e - plot_ts_s + 1

message    = "subday"
pic_name   = "CABLE-HESS_"+message+"_PFT-"+PFT+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"

case_name = (/"CTL-30x3+30yr","FREE_DRAIN", "CTL-30x3+30yr_2019-WetSoil",\
              "CTL-30x3+30yr_2019-DeepRoot"/)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

if PFT .eq. "all" then
  iveg_s = 0
  iveg_e = 20
else if PFT .eq. "tree" then
  iveg_s = 1
  iveg_e = 4
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
  path(case_num) = path_mum+case_name(case_num)+"/outputs-subday"
end do

; ___________________________ TREE Mask ______________________________

filename1 = path(0)+"/cable_out_2019_SE_Aus.nc"
f1 = addfile (filename1,"r")
PFT_mask = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
default_val = f1->Qle@_FillValue

PFT_mask_multi = new((/plot_ts,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, plot_ts-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do

delete(filename1)
delete(f1)

Qle      = new((/case_sum,plot_ts/),float,default_val)
Qh       = new((/case_sum,plot_ts/),float,default_val)
TVeg     = new((/case_sum,plot_ts/),float,default_val)
deltaT   = new((/case_sum,plot_ts/),float,default_val) ; *3
Fwsoil   = new((/case_sum,plot_ts/),float,default_val)

Qle      = 0.
Qh       = 0.
TVeg     = 0.
deltaT   = 0.
Fwsoil   = 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  do year = year_s,year_e
    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc" ; outputs/
    f1 = addfile (filename1,"r")
    ; tmp = f1->Qle(plot_ts_s:plot_ts_e,{A:B},{C:D})/(f1->Qle(plot_ts_s:plot_ts_e,{A:B},{C:D})+f1->Qh(plot_ts_s:plot_ts_e,{A:B},{C:D}))
    ; Qle(case_num,:)      = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1 .and. tmp .gt. 0, tmp, default_val) ,(/1,2/))
    Qle(case_num,:)      = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qle(plot_ts_s:plot_ts_e,{A:B},{C:D}),    default_val) ,(/1,2/))
    Qh(case_num,:)       = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Qh(plot_ts_s:plot_ts_e,{A:B},{C:D}),    default_val) ,(/1,2/))
    TVeg(case_num,:)     = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->TVeg(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/))*3600*24
    deltaT(case_num,:)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->VegT(plot_ts_s:plot_ts_e,{A:B},{C:D}) - f1->Tair(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/))
    ; deltaT(case_num*3,:)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->VegT(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/)) - 273.15
    ; deltaT(case_num*3+1,:) = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Tair(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/)) - 273.15
    ; deltaT(case_num*3+2,:) = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->VegT(plot_ts_s:plot_ts_e,{A:B},{C:D})-f1->Tair(plot_ts_s:plot_ts_e,{A:B},{C:D}),  default_val) ,(/1,2/))
    Fwsoil(case_num,:)   = dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->Fwsoil(plot_ts_s:plot_ts_e,{A:B},{C:D}),default_val) ,(/1,2/))
  end do ; year
  delete(filename1)
  delete(f1)
end do ; case_num
print("prepare to plot")
delete(PFT_mask_multi)

;_____________________________ PLOT ___________________________________
pic = "./plots/"+pic_name
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

res@xyLineThicknesses   = (/  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0,1.0/)          ; make second line thicker
res@xyDashPatterns      = (/0, 0, 0, 0, 0, 0, 0, 0, 0/)
res@xyLineColors        = (/"black","red","blue","green","purple","blue","chartreuse","orange","cyan","yellow"/)          ; change line color
res@xyExplicitLegendLabels = case_name
;res@trYMinF = 0

pres                    = True
res@gsnDraw             = False
res@gsnFrame            = False
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
;pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01


res@tmXBMode  = "Explicit"
x             = ispan(1,plot_ts,1)
res@tmXBValues= ispan(1,plot_ts,8)
res@tmXBLabels= tostring(ispan(11,25,1))
res@tiXAxisString  = "day"
plots = new(5,graphic)
; res@trYMaxF = 0.5
; res@trYMinF = 0.
res@pmLegendDisplayMode = "Never"
res@tiYAxisString  = "Qle (W/m2)"
plots(0) = gsn_csm_xy(wks,x,Qle,res)

res@tiYAxisString  = "Qh (W/m2)"
plots(1) = gsn_csm_xy(wks,x,Qh,res)

res@tiYAxisString  = "deltaT (C)"
plots(2) = gsn_csm_xy(wks,x,deltaT,res)

res@tiYAxisString  = "TVeg (mm/d)"
plots(3) = gsn_csm_xy(wks,x,TVeg,res)

res@pmLegendDisplayMode = "Always"
res@tiYAxisString  = "Fwsoil (-)"
plots(4) = gsn_csm_xy(wks,x,Fwsoil,res)


gsn_panel(wks,(/plots/),(/3,2/),pres)
; txres                   = True
; txres@txPerimOn         = False
; txres@txFont            = 25
; txres@txFontHeightF     = 0.03
;
; amres                   = True
; amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
; amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
; amres@amJust            = "TopLeft"
;
; txid1   = gsn_create_text(wks, "(a)", txres)
; annoid1 = gsn_add_annotation(plot, txid1, amres)
end
