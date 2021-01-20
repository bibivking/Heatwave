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
year_e     = 2019
year_sum   = year_e-year_s+1

; ______________________________ Path _________________________________
path_mum       = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_name      = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr"/) ;,"CTL-30x3+30yr-DeepRoot/outputs"/)
pic_message    = "CTL_vs_FREE_DRAIN"
pic_name       = "HESS_EVAPFBL_"+pic_message+"_"+time_scale+"_"+scale+"_"+year_s+"-"+year_e
Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Depth          = (/0., 1.1, 5.1, 15.7, 43.85, 118.55, 316.4, 460./)

order_name     = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

; "tree" PFT
iveg_s         = 1
iveg_e         = 4

if time_scale .eq. "monthly" then
  print("monthly")
  time_tot = year_sum*12
else if time_scale .eq. "daily" then
  print("daily")
  time_tot = 0
  do year = year_s,year_e
    if mod(year,4) .eq. 0 then
      time_tot = time_tot + 366
    else
      time_tot = time_tot + 365
    end if
  end do
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
filename1   = path(0)+"/cable_out_1990_SE_Aus.nc"
f1          = addfile (filename1,"r")
PFT_mask    = where((f1->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
default_val = f1->Qle@_FillValue

delete(filename1)
delete(f1)

; ___________________________ READ DATA ______________________________
EVAPFBL = new((/case_sum,2,8,time_tot/),float,default_val)
EVAPFBL= 0.

TVeg   = new((/case_sum,2,time_tot/),float,default_val)
TVeg   = 0.

printVarSummary(EVAPFBL)

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

    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
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

    TVeg(case_num, 0, i:i+ts-1)= dim_avg_n_Wrap(f1->TVeg(:,{A:B},{C:D})*1., (/1,2/))
    TVeg(case_num, 1, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->TVeg(:,{A:B},{C:D}),default_val), (/1,2/))
    ; top
    EVAPFBL(case_num, 0, 0, i:i+ts-1)= dim_avg_n_Wrap(f1->EVAPFBL(:,0,{A:B},{C:D})*1., (/1,2/))
    EVAPFBL(case_num, 1, 0, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->EVAPFBL(:,0,{A:B},{C:D}),default_val), (/1,2/))
    do soil = 0,5
       EVAPFBL(case_num, 0, soil+1, i:i+ts-1)= dim_avg_n_Wrap(f1->EVAPFBL(:,soil,{A:B},{C:D})*1., (/1,2/))
       EVAPFBL(case_num, 1, soil+1, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->EVAPFBL(:,soil,{A:B},{C:D}),default_val), (/1,2/))
    end do
    ; bottom
    EVAPFBL(case_num, 0, 7, i:i+ts-1)= dim_avg_n_Wrap(f1->EVAPFBL(:,5,{A:B},{C:D})*1., (/1,2/))
    EVAPFBL(case_num, 1, 7, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->EVAPFBL(:,5,{A:B},{C:D}),default_val), (/1,2/))
    i = i + ts
    delete(PFT_mask_multi)
  end do ; year
end do ; case_num

do ts = 0, time_tot -1
  ; CTL - PFT all
  EVAPFBL(0, 0, :, ts) = EVAPFBL(0, 0, :, ts)/TVeg(0, 0, ts)
  ; CTL - PFT trees
  EVAPFBL(0, 1, :, ts) = EVAPFBL(0, 1, :, ts)/TVeg(0, 1, ts)
  ; FREE_DRAIN - PFT all
  EVAPFBL(1, 0, :, ts) = EVAPFBL(1, 0, :, ts)/TVeg(1, 0, ts)
  ; FREE_DRAIN - PFT trees
  EVAPFBL(1, 1, :, ts) = EVAPFBL(1, 1, :, ts)/TVeg(1, 1, ts)
end do
print("prepare to plot")


; ___________________________ Regrid ______________________________
EVAPFBL!0     = "case"
EVAPFBL!1     = "PFT"
EVAPFBL!2     = "depth"
EVAPFBL!3     = "time"
EVAPFBL&depth = Depth
EVAPFBL&time  = ispan(1,time_tot,1)
Depth_regrid    = fspan(0,200,41) ;fspan(0,460,93)
EVAPFBL_regrid = linint2(EVAPFBL&time,EVAPFBL&depth,EVAPFBL,False,EVAPFBL&time,Depth_regrid,0)
printVarSummary(EVAPFBL_regrid)

EVAPFBL_regrid!0     = "case"
EVAPFBL_regrid!1     = "PFT"
EVAPFBL_regrid!2     = "depth"
EVAPFBL_regrid!3     = "time"
EVAPFBL_regrid&depth = fspan(0,200,41) ;fspan(0,460,93)
EVAPFBL_regrid&time  = ispan(1,time_tot,1)
printVarSummary(EVAPFBL_regrid)
print(EVAPFBL_regrid&depth)

;_____________________________ PLOT ___________________________________
pic = "./plots/"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")

res                         = True
res@gsnDraw                 = False  ; Don't draw plots
res@gsnFrame                = False  ; ����ҳ
res@gsnMaximize             = True

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

res@cnFillOn                = True                            ;��ɫ
res@cnLinesOn               = False
res@cnLineLabelsOn          = False
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode    = "ExplicitLevels"

res@lbLabelBarOn            = True
res@lbBoxLinesOn            = True                       ;�ر�lbar box ����
res@lbTitleFont             = 25
res@lbLabelFont             = 25
res@lbTitleFontHeightF      = 0.015
res@lbLabelFontHeightF      = 0.015
res@lbTitlePosition         = "Right"
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


if time_scale .eq. "monthly" then
  res@tmXBValues            = ispan(1,time_tot,12*2)
  res@tmXBLabels            = ispan(year_s,year_e,2)
  res@tiXAxisString         = "months"
else if time_scale .eq. "daily" then
  res@tmXBValues            = ispan(1,time_tot,366*2)
  res@tmXBLabels            = ispan(year_s,year_e,2)
end if
end if
res@tmYLValues          = fspan(0,200/5,5) ;fspan(10/5,460/5,10)
res@tmYLLabels          = tostring(toint(fspan(200,0,5))) ;tostring(toint(fspan(450,0,10)))


;**************** ͼ�������� *****************
pres                    = True                                      ; ���� panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.2
pres@lbLabelFontHeightF = 0.2

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.03

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"


; _______________________________ plot soil profile _____________________________
plots              = new(6,graphic)

res@cnLevels       = fspan(0.0,50.,21)
res@cnFillPalette  = "WhiteBlueGreenYellowRed";"cmocean_deep";"MPL_YlGn";"GMT_drywet" ;"WhiteBlueGreenYellowRed"
res@tiYAxisString  = "Depth (cm)"
; res@lbTitleString  = "(m3/m3)"
;printVarSummary(EVAPFBL_regrid(0,0,::-1,:)*1.)
plots(0)           = gsn_csm_contour(wks, EVAPFBL_regrid(0,0,::-1,:)*100.,res)
delete(res@tiYAxisString)
plots(1)           = gsn_csm_contour(wks, EVAPFBL_regrid(0,1,::-1,:)*100.,res)
res@tiYAxisString  = "Depth (cm)"
plots(2)           = gsn_csm_contour(wks, EVAPFBL_regrid(1,0,::-1,:)*100.,res)
delete(res@tiYAxisString)
plots(3)           = gsn_csm_contour(wks, EVAPFBL_regrid(1,1,::-1,:)*100.,res)
delete(res@cnLevels)
res@cnLevels       = fspan(-5,5,21)
res@cnFillPalette  = "MPL_PRGn";"cmocean_delta";"MPL_BrBG" ;"ViBlGrWhYeOrRe"
res@tiXAxisString  = "Year"
res@tiYAxisString  = "Depth (cm)"

tmp                =  EVAPFBL_regrid(1,0,::-1,:)*100. - EVAPFBL_regrid(0,0,::-1,:)*100.
printVarSummary(tmp)
; tmp!0              = "depth"
; tmp!1              = "time"
; tmp&depth          = fspan(460,0,93)
; tmp&time           = ispan(1,time_tot,1)
plots(4)           = gsn_csm_contour(wks, tmp*1., res)
delete(res@tiYAxisString)
tmp                = EVAPFBL_regrid(1,1,::-1,:) - EVAPFBL_regrid(0,1,::-1,:)
plots(5)           = gsn_csm_contour(wks, tmp*1., res)

do counter = 0,5
  txid1   = gsn_create_text(wks, order_name(counter), txres)
  annoid1 = gsn_add_annotation(plots(counter), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/3,2/),pres)

end
