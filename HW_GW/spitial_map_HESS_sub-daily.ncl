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
HW_day_s   = 11
HW_day_e   = 26
year_sum   = year_e-year_s+1
time_tot   = year_sum*365*8
plot_ts_s  = (HW_day_s-1)*8-1 ; the first day's 0h
plot_ts_e  = (HW_day_e-1)*8+7 ; the last day's 24h
plot_ts    = plot_ts_e - plot_ts_s + 1

message    = "subday_map"
pic_name   = "CABLE-HESS_"+message+"_PFT-"+PFT+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"

case_name = (/"CTL-30x3+30yr","FREE_DRAIN", "CTL-30x3+30yr_2019-WetSoil",\
              "CTL-30x3+30yr_2019-DeepRoot"/)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers


order_name  = (/"(a)","(b)","(c)","(d)",\
                "(e)","(f)","(g)","(h)",\
                "(i)","(j)","(k)","(l)",\
                "(m)","(n)","(o)","(p)",\
                "(q)","(r)","(s)","(t)"/)

; 0, 1, 2, 3 , 4 , 5,  6,  7,  8
; 0, 3, 6, 9, 12, 15, 18, 21, 24

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
f1  = addfile (filename1,"r")
lat = f1->y({A:B})
lon = f1->x({C:D})

PFT_mask = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
default_val = f1->Qle@_FillValue

PFT_mask_multi = new((/plot_ts,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, plot_ts-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do

delete(filename1)
delete(f1)

Qle      = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val)
Qh       = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val)
TVeg     = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val)
deltaT   = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val) ; *3
Fwsoil   = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val)

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

    Qle_tmp      = f1->Qle(plot_ts_s:plot_ts_e,{A:B},{C:D})
    Qh_tmp       = f1->Qh(plot_ts_s:plot_ts_e,{A:B},{C:D})
    TVeg_tmp     = f1->TVeg(plot_ts_s:plot_ts_e,{A:B},{C:D})
    TVeg_tmp     = TVeg_tmp*3600*3
    deltaT_tmp   = f1->VegT(plot_ts_s:plot_ts_e,{A:B},{C:D})
    deltaT_tmp   = deltaT_tmp - f1->Tair(plot_ts_s:plot_ts_e,{A:B},{C:D})
    Fwsoil_tmp   = f1->Fwsoil(plot_ts_s:plot_ts_e,{A:B},{C:D})
    printVarSummary(Qle_tmp)
    count = 0
    do day = HW_day_s,HW_day_e
      Qle(case_num,:,:)     = Qle(case_num,:,:)    + Qle_tmp(count*8+5,:,:)/(HW_day_e-HW_day_s+1)
      Qh(case_num,:,:)      = Qh(case_num,:,:)     + Qh_tmp(count*8+5,:,:)/(HW_day_e-HW_day_s+1)
      TVeg(case_num,:,:)    = TVeg(case_num,:,:)   + TVeg_tmp(count*8+5,:,:)/(HW_day_e-HW_day_s+1)
      deltaT(case_num,:,:)  = deltaT(case_num,:,:) + deltaT_tmp(count*8+5,:,:)/(HW_day_e-HW_day_s+1)
      Fwsoil(case_num,:,:)  = Fwsoil(case_num,:,:) + Fwsoil_tmp(count*8+5,:,:)/(HW_day_e-HW_day_s+1)
      ; Qle(case_num,:,:)     = Qle(case_num,:,:)    + dim_avg_n_Wrap(Qle_tmp(count*8+3:count*8+7,:,:),0)/(HW_day_e-HW_day_s+1)
      ; Qh(case_num,:,:)      = Qh(case_num,:,:)     + dim_avg_n_Wrap(Qh_tmp(count*8+3:count*8+7,:,:),0)/(HW_day_e-HW_day_s+1)
      ; TVeg(case_num,:,:)    = TVeg(case_num,:,:)   + dim_avg_n_Wrap(TVeg_tmp(count*8+3:count*8+7,:,:),0)/(HW_day_e-HW_day_s+1)
      ; deltaT(case_num,:,:)  = deltaT(case_num,:,:) + dim_avg_n_Wrap(deltaT_tmp(count*8+3:count*8+7,:,:),0)/(HW_day_e-HW_day_s+1)
      ; Fwsoil(case_num,:,:)  = Fwsoil(case_num,:,:) + dim_avg_n_Wrap(Fwsoil_tmp(count*8+3:count*8+7,:,:),0)/(HW_day_e-HW_day_s+1)
      count = count + 1
    end do
  end do ; year
  delete(filename1)
  delete(f1)
end do ; case_num
print("prepare to plot")
delete(PFT_mask_multi)

Qle!0      = "time"
Qle!1      = "lat"
Qle!2      = "lon"
Qle&lat    = lat
Qle&lon    = lon

Qh!0      = "time"
Qh!1      = "lat"
Qh!2      = "lon"
Qh&lat    = lat
Qh&lon    = lon

TVeg!0      = "time"
TVeg!1      = "lat"
TVeg!2      = "lon"
TVeg&lat    = lat
TVeg&lon    = lon

deltaT!0      = "time"
deltaT!1      = "lat"
deltaT!2      = "lon"
deltaT&lat    = lat
deltaT&lon    = lon

Fwsoil!0      = "time"
Fwsoil!1      = "lat"
Fwsoil!2      = "lon"
Fwsoil&lat    = lat
Fwsoil&lon    = lon

; ================================= Plotting ===================================
pic = "./plots/"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")

res                         = True
res@gsnDraw                 = False  ; Don't draw plots
res@gsnFrame                = False  ; ����ҳ
res@gsnMaximize             = True

res@vpHeightF 	            = 0.4               ; change aspect ratio of plot
res@vpWidthF 	              = 0.8
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

res@cnFillMode             = "RasterFill"            ; Raster Mode
res@cnFillOn               = True                            ;��ɫ
res@cnLinesOn              = False
res@cnLineLabelsOn         = False
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode   = "ExplicitLevels"

res@lbLabelBarOn           = True
res@lbBoxLinesOn           = True                       ;�ر�lbar box ����
res@lbTitleFont            = 25
res@lbLabelFont            = 25
res@lbTitleFontHeightF     = 0.013
res@lbLabelFontHeightF     = 0.013

res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015
; res@mpDataBaseVersion       = "HighRes"

res@tmYLMode                = "Explicit"
res@tmXBMode                = "Explicit"
res@txString                = ""
res@gsnLeftString           = ""
res@gsnRightString          = ""


if scale .eq. "Global" then
  res@trYMinF         = -60
  res@trYMaxF         = 80
  res@mpMinLatF       = -60
  res@mpMaxLatF       = 80
  res@mpCenterLonF    = 160
  res@tmYLMode        = "Explicit"
  res@tmYLValues      =(/60,30,0,-30,-60/)
  res@tmYLLabels      =(/"60~S~o~N~N","30~S~o~N~N","0~S~o~N~","30~S~o~N~S","60~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      =(/240,300,0,60,120,180/)
  res@tmXBLabels      =(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)
else if scale .eq. "SE-AU" then
  res@gsnAddCyclic    = False
  res@trYMinF         = -40.
  res@trYMaxF         = -28.
  res@mpMinLatF       = -40.
  res@mpMaxLatF       = -28.
  res@trXMinF         = 140.
  res@trXMaxF         = 154.
  res@mpMinLonF       = 140.
  res@mpMaxLonF       = 154.
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-40,-35,-30,-25/)
  res@tmYLLabels      = (/"40~S~o~N~N","35~S~o~N~N","30~S~o~N~","25~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/140,145,150/)
  res@tmXBLabels      = (/"140~S~o~N~W","145~S~o~N~W","150~S~o~N~"/)
else if scale .eq. "AU" then
  res@gsnAddCyclic    = False
  res@trYMinF         = -44
  res@trYMaxF         = -10
  res@mpMinLatF       = -44
  res@mpMaxLatF       = -10
  res@trXMinF         = 112
  res@trXMaxF         = 154
  res@mpMinLonF       = 112
  res@mpMaxLonF       = 154
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-40,-35,-30,-25,-20,-15,-10/)
  res@tmYLLabels      = (/"40~S~o~N~N","35~S~o~N~N","30~S~o~N~","25~S~o~N~S","20~S~o~N~N","15~S~o~N~N","10~S~o~N~"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/115,120,125,130,135,140,145,150/)
  res@tmXBLabels      = (/"115~S~o~N~W","120~S~o~N~W","125~S~o~N~","130~S~o~N~W","135~S~o~N~W","140~S~o~N~","145~S~o~N~W","150~S~o~N~"/)
end if
end if
end if

;**************** ͼ�������� *****************
pres                    = True
pres@wkWidth            = 2500
pres@wkHeight           = 2500

pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False
pres@cnLevelSelectionMode = "ExplicitLevels"
pres@lbBoxLinesOn       = True
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.05

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

plots = new(case_sum*5,graphic)

; 1: Qle
res@cnLevels      = fspan(0,400,21)
res@cnFillPalette = "WhiteBlueGreenYellowRed";"MPL_RdBu";"vegetation_ClarkU";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
plots(0) = gsn_csm_contour_map(wks,Qle(0,:,:),res)
delete(res@cnLevels)

res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"ViBlGrWhYeOrRe"
res@cnLevels      = fspan(-150,150,21)
tmp  = Qle(1,:,:)
tmp  = tmp - Qle(0,:,:)
plots(1) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = Qle(2,:,:)
tmp  = tmp - Qle(0,:,:)
plots(2) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = Qle(3,:,:)
tmp  = tmp - Qle(0,:,:)
plots(3) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

; do cnt = 0,3
;   txid1   = gsn_create_text(wks, order_name(cnt), txres)
;   annoid1 = gsn_add_annotation(plots(cnt), txid1, amres)
; end do

; gsn_panel(wks,(/plots/),(/1,4/),pres)
delete(res@cnLevels)

count = 4
; 2: Qh
res@cnLevels      = fspan(0,400,21)
res@cnFillPalette = "WhiteBlueGreenYellowRed";"NCV_blue_red";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
plots(0+count) = gsn_csm_contour_map(wks,Qh(0,:,:),res)
delete(res@cnLevels)

res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"ViBlGrWhYeOrRe"
res@cnLevels      = fspan(-150,150,31)
tmp  = Qh(1,:,:)
tmp  = tmp - Qh(0,:,:)
plots(1+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = Qh(2,:,:)
tmp  = tmp - Qh(0,:,:)
plots(2+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = Qh(3,:,:)
tmp  = tmp - Qh(0,:,:)
plots(3+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

; do cnt = 0,3
;   txid1   = gsn_create_text(wks, order_name(cnt), txres)
;   annoid1 = gsn_add_annotation(plots(cnt+count), txid1, amres)
; end do

; gsn_panel(wks,(/plots/),(/1,4/),pres)
delete(res@cnLevels)


count = 8
; 3: TVeg
res@cnLevels      = fspan(0,2,21)
res@cnFillPalette = "WhiteBlueGreenYellowRed";"MPL_RdBu";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
plots(0+count) = gsn_csm_contour_map(wks,TVeg(0,:,:),res)
delete(res@cnLevels)

res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"ViBlGrWhYeOrRe"
res@cnLevels      = fspan(-1.,1.,21)
tmp  = TVeg(1,:,:)
tmp  = tmp - TVeg(0,:,:)
plots(1+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = TVeg(2,:,:)
tmp  = tmp - TVeg(0,:,:)
plots(2+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = TVeg(3,:,:)
tmp  = tmp - TVeg(0,:,:)
plots(3+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

; do cnt = 0,3
;   txid1   = gsn_create_text(wks, order_name(cnt), txres)
;   annoid1 = gsn_add_annotation(plots(cnt+count), txid1, amres)
; end do

; gsn_panel(wks,(/plots/),(/1,4/),pres)
delete(res@cnLevels)

count = 12
; 4: deltaT
res@cnLevels      = fspan(0,15,16)
res@cnFillPalette = "WhiteBlueGreenYellowRed";"NCV_blue_red";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
plots(0+count) = gsn_csm_contour_map(wks,deltaT(0,:,:),res)
delete(res@cnLevels)


res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"ViBlGrWhYeOrRe"
res@cnLevels      = fspan(-5.,5.,21)
tmp  = deltaT(1,:,:)
tmp  = tmp - deltaT(0,:,:)
plots(1+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = deltaT(2,:,:)
tmp  = tmp - deltaT(0,:,:)
plots(2+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = deltaT(3,:,:)
tmp  = tmp - deltaT(0,:,:)
plots(3+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

; do cnt = 0,3
;   txid1   = gsn_create_text(wks, order_name(cnt), txres)
;   annoid1 = gsn_add_annotation(plots(cnt+count), txid1, amres)
; end do

; gsn_panel(wks,(/plots/),(/1,4/),pres)
delete(res@cnLevels)

count = 16
; 5: Fwsoil
res@cnLevels      = fspan(0,1.,21)
res@cnFillPalette = "WhiteBlueGreenYellowRed";"MPL_RdBu";"vegetation_ClarkU";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
plots(0+count) = gsn_csm_contour_map(wks,Fwsoil(0,:,:),res)
delete(res@cnLevels)


res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"ViBlGrWhYeOrRe"
res@cnLevels      = fspan(-0.5,0.5,21)
tmp  = Fwsoil(1,:,:)
tmp  = tmp - Fwsoil(0,:,:)
plots(1+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = Fwsoil(2,:,:)
tmp  = tmp - Fwsoil(0,:,:)
plots(2+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

tmp  = Fwsoil(3,:,:)
tmp  = tmp - Fwsoil(0,:,:)
plots(3+count) = gsn_csm_contour_map(wks,tmp,res)
delete(tmp)

do cnt = 0,19
  txid1   = gsn_create_text(wks, order_name(cnt), txres)
  annoid1 = gsn_add_annotation(plots(cnt), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/5,4/),pres)

end
