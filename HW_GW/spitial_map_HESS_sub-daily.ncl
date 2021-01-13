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
PFT        = "all" ; "tree"
scale      = "SE-AU"
year_s     = 2019
year_e     = 2019

HW_day_s   = (/ 1 /) ;15, 22/)
HW_day_e   = (/ 31 /) ;19, 27/)
event_sum  = dimsizes(HW_day_s)

day_tot    = 0
do event_num = 0, event_sum-1
  day_tot = day_tot + (HW_day_e(event_num)-HW_day_s(event_num)+1)
end do

message    = "subday_map"
pic_name   = "CABLE-HESS_"+message+"_PFT-"+PFT+"_"+scale+"_"+year_s+"-"+year_e

path_mum  = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"

case_name = (/"CTL-30x3+30yr","FREE_DRAIN", "CTL-30x3+30yr_2019-WetSoil",\
              "CTL-30x3+30yr_2019-DeepRoot"/)

order_name  = (/"(a)","(b)","(c)","(d)",\
                "(e)","(f)","(g)","(h)",\
                "(i)","(j)","(k)","(l)",\
                "(m)","(n)","(o)","(p)",\
                "(q)","(r)","(s)","(t)"/)

; 0, 1, 2, 3 , 4 , 5,  6,  7
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

delete(filename1)
delete(f1)

Qle      = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val)
Qh       = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val)
; TVeg     = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val)
deltaT   = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val) ; *3
Fwsoil   = new((/case_sum,dimsizes(lat),dimsizes(lon)/),float,default_val)

Qle      = 0.
Qh       = 0.
; TVeg     = 0.
deltaT   = 0.
Fwsoil   = 0.

do case_num = 0, case_sum -1
  print("case_num = " +case_num)
  do year = year_s,year_e
    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc" ; outputs/
    f1 = addfile (filename1,"r")

    Qle_tmp      = f1->Qle(0:60*8,{A:B},{C:D})
    Qh_tmp       = f1->Qh(0:60*8,{A:B},{C:D})
    deltaT_tmp   = f1->VegT(0:60*8,{A:B},{C:D})
    deltaT_tmp   = deltaT_tmp - f1->Tair(0:60*8,{A:B},{C:D})
    Fwsoil_tmp   = f1->Fwsoil(0:60*8,{A:B},{C:D})
    printVarSummary(Qle_tmp)
    delete(f1)
    delete(filename1)

    do event_num = 0,event_sum-1
      do day = HW_day_s(event_num),HW_day_e(event_num)
        ; ts = (HW_day_s(event_num)-1)*8+4
        ts = (day-1)*8+6;6-6pm ;7 - 9pm
        Qle(case_num,:,:)    = Qle(case_num,:,:) + where(PFT_mask .eq. 1, Qle_tmp(ts,:,:) - Qle_tmp(ts+1,:,:), default_val)/day_tot
        Qh(case_num,:,:)     = Qh(case_num,:,:)  + where(PFT_mask .eq. 1,Qh_tmp(ts,:,:) - Qh_tmp(ts+1,:,:), default_val)/day_tot
        deltaT(case_num,:,:) = deltaT(case_num,:,:) + where(PFT_mask .eq. 1,deltaT_tmp(ts,:,:) - deltaT_tmp(ts+1,:,:), default_val)/day_tot
        Fwsoil(case_num,:,:) = Fwsoil(case_num,:,:) + where(PFT_mask .eq. 1,Fwsoil_tmp(ts,:,:) - Fwsoil_tmp(ts+1,:,:), default_val)/day_tot
      end do
    end do
  end do ; year
end do ; case_num
delete(PFT_mask)

print("prepare to plot")
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

; TVeg!0      = "time"
; TVeg!1      = "lat"
; TVeg!2      = "lon"
; TVeg&lat    = lat
; TVeg&lon    = lon

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
res@gsnMaximize             = False

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
; res@lbOrientation          = "vertical"
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
  res@tmYLLabels      = (/"40~S~o~N~S","35~S~o~N~S","30~S~o~N~S","25~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/140,145,150/)
  res@tmXBLabels      = (/"140~S~o~N~E","145~S~o~N~E","150~S~o~N~E"/)
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
  res@tmYLLabels      = (/"40~S~o~N~S","35~S~o~N~S","30~S~o~N~S","25~S~o~N~S","20~S~o~N~S","15~S~o~N~S","10~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/115,120,125,130,135,140,145,150/)
  res@tmXBLabels      = (/"115~S~o~N~E","120~S~o~N~E","125~S~o~N~E","130~S~o~N~E","135~S~o~N~E","140~S~o~N~E","145~S~o~N~E","150~S~o~N~E"/)
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
pres@gsnPanelYWhiteSpacePercent = 0.5
pres@gsnPanelXWhiteSpacePercent = 0.5

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.03

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

Qle_CTL_lvl    = (/25., 50, 100, 150, 200, 250, 300, 350, 400/)
Fw_CTL_lvl     = fspan(0.1,0.9,9)
Qh_CTL_lvl     = Qle_CTL_lvl
deltaT_CTL_lvl = fspan(1,12,12)

; diff_color = "BlAqGrWh2YeOrReVi22"; "ViBlGrWhYeOrRe";
Qle_diff_lvl    = (/-200, -150, -100, -50, -25., 0, 25., 50, 100, 150, 200/) ;-40, -30, -20, -10;10, 20, 30, 40,
Fw_diff_lvl     = fspan(-0.6,0.6,13)
Qh_diff_lvl     = Qle_diff_lvl
deltaT_diff_lvl = (/-5., -4., -3., -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4., 5./)

plots = new(case_sum*4,graphic)
do case_num = 0,case_sum-1
  ; if case_num .eq. 0 then
    res@lbLabelBarOn  = True

    ; deltaT
    res@cnLevels      = deltaT_diff_lvl ;deltaT_CTL_lvl
    res@cnFillPalette = "BlAqGrWh2YeOrReVi22"
    ; res@cnFillPalette = "precip3_16lev";"WhiteBlueGreenYellowRed";"vegetation_ClarkU";"NCV_blue_red";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
    ; plots(case_num*4) = gsn_csm_contour_map(wks,deltaT(0,:,:),res)
    plots(case_num*4) = gsn_csm_contour_map(wks,deltaT(case_num,:,:),res)
    delete(res@cnLevels)

    ; Qh
    res@cnLevels      = Qh_diff_lvl ;Qh_CTL_lvl
    ; plots(case_num*4+1) = gsn_csm_contour_map(wks,Qh(0,:,:),res)
    plots(case_num*4+1) = gsn_csm_contour_map(wks,Qh(case_num,:,:),res)
    delete(res@cnLevels)
    delete(res@cnFillPalette)

    ; Qle
    delete(res@cnFillPalette)
    cmap = read_colormap_file("BlAqGrWh2YeOrReVi22");("NEO_div_vegetation_a")("ViBlGrWhYeOrRe");
    res@cnFillPalette = cmap(::-1,:)

    res@cnLevels      = Qle_diff_lvl ;Qle_CTL_lvl
    ;cmap = read_colormap_file("WhiteBlueGreenYellowRed");("NEO_div_vegetation_a")("ViBlGrWhYeOrRe");
    ; res@cnFillPalette = "spread_15lev";cmap(::-1,:)
    ; plots(case_num*4+2) = gsn_csm_contour_map(wks,Qle(0,:,:),res)
    plots(case_num*4+2) = gsn_csm_contour_map(wks,Qle(case_num,:,:),res)
    delete(res@cnLevels)

    ; Fwsoil
    res@cnLevels      = fspan(-0.1,0.1,21) ;Fw_CTL_lvl
    ; plots(case_num*4+3) = gsn_csm_contour_map(wks,Fwsoil(0,:,:),res)
    plots(case_num*4+3) = gsn_csm_contour_map(wks,Fwsoil(case_num,:,:),res)
    delete(res@cnLevels)
    delete(res@cnFillPalette)

  ; else
  ;   if case_num .lt.case_sum-1 then
  ;     res@lbLabelBarOn  = False
  ;   else
  ;     res@lbLabelBarOn  = True
  ;   end if
  ;   delete(res@cnFillPalette)
  ;   res@cnLevels      = deltaT_diff_lvl
  ;   res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"BlAqGrWh2YeOrReVi22"
  ;   tmp  = deltaT(case_num,:,:)
  ;   tmp  = tmp - deltaT(0,:,:)
  ;   plots(case_num*4) = gsn_csm_contour_map(wks,tmp,res)
  ;   delete(tmp)
  ;   delete(res@cnLevels)
  ;
  ;   res@cnLevels      = Qh_diff_lvl
  ;   tmp  = Qh(case_num,:,:)
  ;   tmp  = tmp - Qh(0,:,:)
  ;   plots(case_num*4+1) = gsn_csm_contour_map(wks,tmp,res)
  ;   delete(tmp)
  ;   delete(res@cnLevels)
  ;
  ;   ; delete(cmap)
  ;   delete(res@cnFillPalette)
  ;   cmap = read_colormap_file("BlAqGrWh2YeOrReVi22");("NEO_div_vegetation_a")("ViBlGrWhYeOrRe");
  ;   res@cnFillPalette = cmap(::-1,:)
  ;
  ;   res@cnLevels      = Qle_diff_lvl
  ;   tmp  = Qle(case_num,:,:)
  ;   tmp  = tmp - Qle(0,:,:)
  ;   plots(case_num*4+2) = gsn_csm_contour_map(wks,tmp,res)
  ;   delete(tmp)
  ;   delete(res@cnLevels)
  ;
  ;   res@cnLevels      = Fw_diff_lvl
  ;   tmp  = Fwsoil(case_num,:,:)
  ;   tmp  = tmp - Fwsoil(0,:,:)
  ;   plots(case_num*4+3) = gsn_csm_contour_map(wks,tmp,res)
  ;   delete(tmp)
  ;   delete(res@cnLevels)
  ;
  ; end if
end do

do cnt = 0, 15
  txid1   = gsn_create_text(wks, order_name(cnt), txres)
  annoid1 = gsn_add_annotation(plots(cnt), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/4,4/),pres)

end


; ; 3: TVeg
; res@cnLevels      = fspan(0,2,21)
; res@cnFillPalette = "WhiteBlueGreenYellowRed";"MPL_RdBu";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
; plots(0+count) = gsn_csm_contour_map(wks,TVeg(0,:,:),res)
; delete(res@cnLevels)
;
; res@cnFillPalette = "BlAqGrWh2YeOrReVi22";"ViBlGrWhYeOrRe"
; res@cnLevels      = fspan(-1.,1.,21)
; tmp  = TVeg(1,:,:)
; tmp  = tmp - TVeg(0,:,:)
; plots(1+count) = gsn_csm_contour_map(wks,tmp,res)
; delete(tmp)
;
; tmp  = TVeg(2,:,:)
; tmp  = tmp - TVeg(0,:,:)
; plots(2+count) = gsn_csm_contour_map(wks,tmp,res)
; delete(tmp)
;
; tmp  = TVeg(3,:,:)
; tmp  = tmp - TVeg(0,:,:)
; plots(3+count) = gsn_csm_contour_map(wks,tmp,res)
; delete(tmp)
; delete(res@cnLevels)
;
; count = 12
