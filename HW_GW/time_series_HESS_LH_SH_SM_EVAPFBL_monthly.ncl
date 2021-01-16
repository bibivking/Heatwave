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
year_s       = 2000
year_e       = 2019
scale        = "SE-AU"

time_scale   = "daily" ;"monthly" ; "daily"
LH_SH_scale  = "seasonal"
EVAPFBL_plot = "bar" ;"line" ;"profile" ;
EVAPFBL_scale= "annual"
               ; options for line plot: "seasonal"; "monthly"
               ;         for bar  plot: "annual"

year_sum     = year_e-year_s+1
if EVAPFBL_plot .eq. "line" .or. EVAPFBL_plot .eq. "bar" then
  EVAPFBL_deseasonal = False
end if

; ______________________________ Path _________________________________
path_mum       = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_name      = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr",\
                   "CTL-30x3+30yr-DeepRoot/outputs-daily"/);,"FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
pic_message    = "FD_DR_CTL_diff_trees"

pic_name       = "HESS_LH_SH_SM_EVAPFBL-"+EVAPFBL_plot+"-"+EVAPFBL_scale+"_"+pic_message+"_"+time_scale+"_"+scale+"_"+year_s+"-"+year_e
Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Depth          = (/1.1, 5.1, 15.7, 43.85, 118.55, 316.4, 2355.36/)
                 ; avg(GWdz) over SE Aus = 37.9072 m
                 ; -> last depth=460+(37.9072*10)/2 = 2355.36 (cm)

Depth1         = (/0., 1.1, 5.1, 15.7, 43.85, 118.55, 316.4, 460./) ; for EVAPFBL
layer_depth    = (/"lvl 1","lvl 2","lvl 3","lvl 4","lvl 5","lvl 6"/);,\
                  ; "L1","L2","L3","L4","L5","L6"/)
                 ;(/"0~2.2cm","2.2~8cm","8~23.4cm","23.4~64.3cm","64.3~172.8cm","172.8~460cm"/)

order_name     = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

iveg_tree_s    = 1
iveg_tree_e    = 4
iveg_nontree_s = 5
iveg_nontree_e = 17

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
f1          = addfile(filename1,"r")
PFT_mask    = where((f1->iveg({A:B},{C:D}) .ge. iveg_tree_s) .and. (f1->iveg({A:B},{C:D}) .le. iveg_tree_e), 1, f1->iveg({A:B},{C:D}))
PFT_mask    = where((PFT_mask .ge. iveg_nontree_s) .and. (PFT_mask .le. iveg_nontree_e), 2, PFT_mask)
default_val = f1->Qle@_FillValue
PFT_mask_soil = new((/6,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, 5
  PFT_mask_soil(z,:,:) = PFT_mask
end do

delete(filename1)
delete(f1)
; ___________________________ READ DATA ______________________________
SoilMoist = new((/case_sum,2,7,time_tot/),float,default_val)
Qle       = new((/case_sum,2,year_sum*12/),float,default_val)
Qh        = new((/case_sum,2,year_sum*12/),float,default_val)
TVeg      = new((/case_sum,2,year_sum*12/),float,default_val)

if EVAPFBL_plot .eq. "profile" then
  EVAPFBL   = new((/case_sum,2,8,time_tot/),float,default_val)
else if EVAPFBL_plot .eq. "line" .or. EVAPFBL_plot .eq. "bar" then
  EVAPFBL   = new((/case_sum,2,6,year_sum*12/),float,default_val)
end if
end if

SoilMoist = 0.
Qle       = 0.
Qh        = 0.
TVeg      = 0.
EVAPFBL   = 0.

printVarSummary(SoilMoist)

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

    Qle_tmp  = f1->Qle(:,{A:B},{C:D})
    Qh_tmp   = f1->Qh(:,{A:B},{C:D})
    TVeg_tmp = f1->TVeg(:,{A:B},{C:D})
    EVAPFBL_tmp  = f1->EVAPFBL(:,:,{A:B},{C:D})

    do month_num = 0,11
      Qle(case_num, 0, j)  = avg(Qle_tmp(dom(month_num):dom(month_num+1)-1,:,:))
      Qh(case_num, 0, j)   = avg(Qh_tmp(dom(month_num):dom(month_num+1)-1,:,:))
      TVeg(case_num, 0, j) = avg(dim_sum_n_Wrap(TVeg_tmp(dom(month_num):dom(month_num+1)-1,:,:), 0))

      Qle(case_num, 1, j)  = avg(where(PFT_mask .eq. 1, \
                                       dim_avg_n_Wrap(Qle_tmp(dom(month_num):dom(month_num+1)-1,:,:), 0),\
                                       default_val))
      Qh(case_num, 1, j)   = avg(where(PFT_mask .eq. 1, \
                                       dim_avg_n_Wrap(Qh_tmp(dom(month_num):dom(month_num+1)-1,:,:), 0),\
                                       default_val))
      TVeg(case_num, 1, j) = avg(where(PFT_mask .eq. 1, \
                                       dim_sum_n_Wrap(TVeg_tmp(dom(month_num):dom(month_num+1)-1,:,:), 0),\
                                       default_val))
      if EVAPFBL_plot .eq. "line" .or. EVAPFBL_plot .eq. "bar" then
         ; print("j="+j)
         EVAPFBL(case_num, 0, :, j)  = dim_avg_n( \
                                       dim_sum_n_Wrap( EVAPFBL_tmp(dom(month_num):dom(month_num+1)-1,:,:,:), 0), \
                                       (/1,2/))
         EVAPFBL(case_num, 1, :, j)  = dim_avg_n(where(PFT_mask_soil .eq. 1, \
                                       dim_sum_n_Wrap(EVAPFBL_tmp(dom(month_num):dom(month_num+1)-1,:,:,:), 0),\
                                       default_val), (/1,2/))
      end if

      j = j + 1
    end do ; month_num
    delete(Qle_tmp)
    delete(Qh_tmp)
    delete(TVeg_tmp)
    delete(EVAPFBL_tmp)

    PFT_mask_multi = new((/ts,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
    do z = 0, ts-1
      PFT_mask_multi(z,:,:) = PFT_mask
    end do

    SoilMoist_tmp = f1->SoilMoist(:,:,{A:B},{C:D})
    do soil = 0,5
       SoilMoist(case_num, 0, soil, i:i+ts-1)= dim_avg_n_Wrap(SoilMoist_tmp(:,soil,:,:)*1., (/1,2/))
       SoilMoist(case_num, 1, soil, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, SoilMoist_tmp(:,soil,:,:),default_val), (/1,2/))
    end do
    delete(SoilMoist_tmp)
    ; aquifer
    SoilMoist(case_num, 0, 6, i:i+ts-1)= dim_avg_n_Wrap(f1->GWMoist(:,{A:B},{C:D})*1., (/1,2/))
    SoilMoist(case_num, 1, 6, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, f1->GWMoist(:,{A:B},{C:D}),default_val), (/1,2/))

    if EVAPFBL_plot .eq. "profile" then
      EVAPFBL_tmp = f1->EVAPFBL(:,:,{A:B},{C:D})
      ; top
      EVAPFBL(case_num, 0, 0, i:i+ts-1)= dim_avg_n_Wrap(EVAPFBL_tmp(:,0,:,:)*1., (/1,2/))
      EVAPFBL(case_num, 1, 0, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, EVAPFBL_tmp(:,0,:,:),default_val), (/1,2/))
      do soil = 0,5
         EVAPFBL(case_num, 0, soil+1, i:i+ts-1)= dim_avg_n_Wrap(EVAPFBL_tmp(:,soil,:,:)*1., (/1,2/))
         EVAPFBL(case_num, 1, soil+1, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, EVAPFBL_tmp(:,soil,:,:),default_val), (/1,2/))
      end do
      ; bottom
      EVAPFBL(case_num, 0, 7, i:i+ts-1)= dim_avg_n_Wrap(EVAPFBL_tmp(:,5,:,:)*1., (/1,2/))
      EVAPFBL(case_num, 1, 7, i:i+ts-1)= dim_avg_n_Wrap(where(PFT_mask_multi .eq. 1, EVAPFBL_tmp(:,5,:,:),default_val), (/1,2/))
      delete(EVAPFBL_tmp)
    end if

    i = i + ts
    delete(PFT_mask_multi)
  end do ; year
end do ; case_num
TVeg    = TVeg*24.*3600.
EVAPFBL = EVAPFBL*24.*3600.

TMP     = dim_sum_n_Wrap(EVAPFBL,2)
TMP     = TMP - TVeg
print(dim_sum_n_Wrap(TMP,2))
delete(TMP)

if EVAPFBL_deseasonal then
  EVAPFBL = rmMonAnnCycLLLT(EVAPFBL)
end if
print("prepare to plot")

; ___________________________ Regrid ______________________________
SoilMoist!0     = "case"
SoilMoist!1     = "PFT"
SoilMoist!2     = "depth"
SoilMoist!3     = "time"
SoilMoist&depth = Depth
SoilMoist&time  = ispan(1,time_tot,1)
Depth_regrid    = fspan(0,460,93)
SoilMoist_regrid = linint2(SoilMoist&time,SoilMoist&depth,SoilMoist,False,SoilMoist&time,Depth_regrid,0)
printVarSummary(SoilMoist_regrid)

SoilMoist_regrid!0     = "case"
SoilMoist_regrid!1     = "PFT"
SoilMoist_regrid!2     = "depth"
SoilMoist_regrid!3     = "time"
SoilMoist_regrid&depth = fspan(0,460,93)
SoilMoist_regrid&time  = ispan(1,time_tot,1)
printVarSummary(SoilMoist_regrid)
; print(SoilMoist_regrid&depth)

if EVAPFBL_plot .eq. "profile" then
  printVarSummary(EVAPFBL)
  EVAPFBL!0     = "case"
  EVAPFBL!1     = "PFT"
  EVAPFBL!2     = "depth"
  EVAPFBL!3     = "time"
  EVAPFBL&depth = Depth1
  EVAPFBL&time  = ispan(1,time_tot,1)
  Depth_regrid  = fspan(0,460,93)
  EVAPFBL_regrid = linint2(EVAPFBL&time,EVAPFBL&depth,EVAPFBL,False,EVAPFBL&time,Depth_regrid,0)
  printVarSummary(EVAPFBL_regrid)

  EVAPFBL_regrid!0     = "case"
  EVAPFBL_regrid!1     = "PFT"
  EVAPFBL_regrid!2     = "depth"
  EVAPFBL_regrid!3     = "time"
  EVAPFBL_regrid&depth = fspan(0,460,93)
  EVAPFBL_regrid&time  = ispan(1,time_tot,1)
  printVarSummary(EVAPFBL_regrid)
  print(EVAPFBL_regrid&depth)
end if

; _________________________ LH SH monthly to seasonaly _________________________
if LH_SH_scale .eq. "seasonal" then
  Qle!0 = "case"
  Qle!1 = "PFT"
  Qle!2 = "time"
  Qh!0  = "case"
  Qh!1  = "PFT"
  Qh!2  = "time"
  Qle_season  = new((/case_sum,2,year_sum*4/),float, default_val)
  Qh_season   = new((/case_sum,2,year_sum*4/),float, default_val)
  do case_num = 0, case_sum-1
    do pft = 0, 1
      print(dimsizes(month_to_seasonN(Qle(case_num,pft,:),(/"DJF","MAM","JJA","SON"/))))
      tmp_Qle = month_to_seasonN(Qle(case_num,pft,:),(/"DJF","MAM","JJA","SON"/))
      tmp_Qle!0 = "season"
      tmp_Qle!1 = "year"
      Qle_season(case_num,pft,:) = ndtooned(tmp_Qle(year|:,season|:))
      delete(tmp_Qle)
      tmp_Qh = month_to_seasonN(Qh(case_num,pft,:),(/"DJF","MAM","JJA","SON"/))
      tmp_Qh!0 = "season"
      tmp_Qh!1 = "year"
      Qh_season(case_num,pft,:)  = ndtooned(tmp_Qh(year|:,season|:))
      delete(tmp_Qh)
    end do
  end do
  printVarSummary(Qle_season)
end if

; _________________________ EVAPFBL monthly to seasonaly _________________________
if EVAPFBL_scale .eq. "seasonal" then
  printVarSummary(EVAPFBL)
  EVAPFBL!0 = "case"
  EVAPFBL!1 = "PFT"
  EVAPFBL!2 = "soil"
  EVAPFBL!3 = "time"

  EVAPFBL_tmp = month_to_seasonN(EVAPFBL(time|:, case|:,PFT|:, soil|:),(/"DJF","MAM","JJA","SON"/))
  printVarSummary(EVAPFBL_tmp)
  EVAPFBL_tmp!0 = "season"
  EVAPFBL_tmp!1 = "year"
  EVAPFBL_tmp!2 = "case"
  EVAPFBL_tmp!3 = "PFT"
  EVAPFBL_tmp!4 = "soil"

  EVAPFBL_season  = new((/case_sum,2,6,year_sum*4/),float, default_val)
  do case_num = 0, case_sum-1
    do pft = 0, 1
      do soil = 0, 5
        EVAPFBL_season(case_num,pft,soil,:) = ndtooned(EVAPFBL_tmp(year|:,season|:,case|case_num,PFT|pft,soil|soil))
      end do
    end do
  end do
  delete(EVAPFBL_tmp)
  printVarSummary(EVAPFBL_season)
else if EVAPFBL_scale .eq. "annual" then
  printVarSummary(EVAPFBL)
  EVAPFBL!0 = "case"
  EVAPFBL!1 = "PFT"
  EVAPFBL!2 = "soil"
  EVAPFBL!3 = "time"

  EVAPFBL_tmp = month_to_annual(EVAPFBL(time|:, case|:,PFT|:, soil|:),0)
  printVarSummary(EVAPFBL_tmp)
  EVAPFBL_tmp!0 = "year"
  EVAPFBL_tmp!1 = "case"
  EVAPFBL_tmp!2 = "PFT"
  EVAPFBL_tmp!3 = "soil"

  EVAPFBL_annual  = new((/case_sum,2,2,year_sum/),float, default_val)
  do case_num = 0, case_sum-1
    do pft = 0, 1
      EVAPFBL_annual(case_num,pft,0,:) = dim_sum_n_Wrap(EVAPFBL_tmp(case|case_num,PFT|pft,soil|:,year|:),0)
      EVAPFBL_annual(case_num,pft,1,:) = dim_sum_n_Wrap(EVAPFBL_tmp(case|case_num,PFT|pft,soil|4:5,year|:),0)
    end do
  end do
  delete(EVAPFBL_tmp)
  printVarSummary(EVAPFBL_annual)

  pfts    = 2
  bars    = 2 ;(CTL FD & CTL DR)
  depthes = 2 ;(all layers & deep 2 layers)
  EVAPFBL_bar    = new((/pfts,bars,depthes,year_sum*2/),float, default_val)
  counter = 0

  do year = year_s -year_s, year_e-year_s
    EVAPFBL_bar(:,0,:,counter)   = EVAPFBL_annual(0,:,:,year)
    EVAPFBL_bar(:,0,:,counter+1) = EVAPFBL_annual(1,:,:,year)
    EVAPFBL_bar(:,1,:,counter)   = EVAPFBL_annual(0,:,:,year)
    EVAPFBL_bar(:,1,:,counter+1) = EVAPFBL_annual(2,:,:,year)
    counter = counter + 2
  end do
  delete(EVAPFBL_annual)
end if
end if

;_____________________________ PLOT ___________________________________
pic = "./plots/"+pic_name
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")

res                         = True
res@gsnDraw                 = True  ; Don't draw plots
res@gsnFrame                = True  ; ����ҳ
res@gsnMaximize             = False

res@vpHeightF               = 0.3
res@vpWidthF                = 0.6
res@tiYAxisFontThicknessF   = 2.5
res@tmXBMajorThicknessF     = 2.5
res@tmXBMinorThicknessF     = 2.5
res@tmYLMajorThicknessF     = 2.5
res@tmYLMinorThicknessF     = 2.5
res@tmYLLabelFontThicknessF = 2.5
res@tmBorderThicknessF      = 2.5
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

res1                        = res
res2                        = res

res@cnFillOn                = True                            ;��ɫ
res@cnLinesOn               = False
res@cnLineLabelsOn          = False
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode    = "ExplicitLevels"

res@lbLabelBarOn            = True
res@cnInfoLabelOn           = False
; res@cnInfoLabelString       = ""
res@lbBoxLinesOn            = False                       ;�ر�lbar box ����
res@lbLabelFont             = 25
res@lbOrientation           = "Vertical"
res@lbLabelFontHeightF      = 0.015
res@lbTitleFont             = 25
res@lbTitlePosition         = "Bottom"
res@lbTitleFontHeightF      = 0.015
res@trXMinF                 = 1
res@trXMaxF                 = time_tot+1

res1@lgPerimOn              = False      ; No legend perimeter.
res1@lgJustification        = "TopRight"
res1@pmLegendParallelPosF   = 0.5       ; Move to right
res1@pmLegendOrthogonalPosF = 0.5     ; Move into plot
res1@pmLegendWidthF         = 0.3       ; Decrease width
res1@pmLegendHeightF        = 0.1       ; Decrease height
res1@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines

if time_scale .eq. "monthly" then
  res@tmXBValues            = ispan(1,time_tot+1,12*2)
  res@tmXBLabels            = ispan(year_s,year_e+1,2)
  res@tiXAxisString         = "months"
else if time_scale .eq. "daily" then
  res@tmXBValues            = ispan(1,time_tot,365*2)
  res@tmXBLabels            = ispan(year_s,year_e+1,2)
end if
end if
res@tmYLValues          = fspan(10/5,460/5,10)
res@tmYLLabels          = tostring(toint(fspan(450,0,10)))


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
pres@gsnPanelYWhiteSpacePercent = 0.1
pres@gsnPanelXWhiteSpacePercent = 0.1
pres@gsnPaperOrientation = "landscape" ; This resource only comes into effect if gsnMaximize is set to True
pres@gsnPaperWidth      = 0.5

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.03

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

plots                   = new(6,graphic)

; ________________________________ Plot Qh Qle _________________________________
; res@vpXF                    = 0.03
; res@vpYF                    = 0.7
res1@gsnYRefLine    = 0.              ; reference line
res1@gsnXYBarChart  = True            ; create bar chart
res1@trYMinF        = -18.
res1@trYMaxF        = 18.
res1@pmLegendDisplayMode  = "Never"

if LH_SH_scale .eq. "seasonal" then
  season_tot         = year_sum*4
  res1@tmXBValues  = ispan(1,season_tot+1,4*2)
  x                = ispan(1,season_tot,1)-1./5.
  res1@tmXBLabels  = ispan(year_s,year_e+1,2)
  res1@trXMinF     = 1.
  res1@trXMaxF     = season_tot+1

  ; ; all
  ; res1@tiYAxisString        = "Heat flux (W/m2)"
  ; res1@pmLegendDisplayMode  = "Never"
  ; res1@gsnXYBarChartColors2 = (/"blue","white","white","white"/) ; gsnBelowYRefLineBarColors
  ; res1@xyLineColor          = "blue" ;(/"blue","black","black","black"/)
  ; plot_Qle                  = gsn_csm_xy(wks, x, Qle_season(1,0,:)-Qle_season(0,0,:),res1)
  ; delete(res1@tiYAxisString)

  ; res1@gsnXYBarChartColors2 = (/"red","white","white","white"/)
  ; res1@xyLineColor          = "red" ;(/"red","black","black","black"/)
  ; plot_Qh                   = gsn_csm_xy(wks, x, Qh_season(1,0,:)-Qh_season(0,0,:),res1)
  ; overlay(plot_Qh, plot_Qle)
  ;
  ; plots(0) = plot_Qh ;plot_Qh

  ; trees
  res1@gsnXYBarChartColors2 = (/"blue","white","white","white"/) ; gsnBelowYRefLineBarColors
  res1@xyLineColor          = "blue" ;(/"blue","black","black","black"/)
  ; FD - CTL
  res1@tiYAxisString        = "~F5~D"+"~F25~LH, "+"~F5~D"+"~F25~SH (W m~S~-2~N~)"
  plot_Qle1                 = gsn_csm_xy(wks, x, Qle_season(1,1,:)-Qle_season(0,1,:),res1)
  delete(res1@tiYAxisString)

  ; DR - CTL
  plot_Qle2                 = gsn_csm_xy(wks, x, Qle_season(2,1,:)-Qle_season(0,1,:),res1)
  ; ; FD+DR - FD
  ; plot_Qle3                 = gsn_csm_xy(wks, x, Qle_season(3,1,:)-Qle_season(1,1,:),res1)

  res1@gsnXYBarChartColors2 = (/"red","white","white","white"/)
  res1@xyLineColor          = "red" ;(/"red","black","black","black"/)
  ; FD - CTL
  res1@tiYAxisString       = "~F5~D"+"~F25~LH, "+"~F5~D"+"~F25~SH (W m~S~-2~N~)"
  plot_Qh1                 = gsn_csm_xy(wks, x, Qh_season(1,1,:)-Qh_season(0,1,:),res1)
  delete(res1@tiYAxisString)
  ; DR - CTL
  plot_Qh2                 = gsn_csm_xy(wks, x, Qh_season(2,1,:)-Qh_season(0,1,:),res1)
  ; ; FD+DR - FD
  ; plot_Qh3                 = gsn_csm_xy(wks, x, Qh_season(3,1,:)-Qh_season(1,1,:),res1)
  overlay(plot_Qh1, plot_Qle1)
  overlay(plot_Qh2, plot_Qle2)
  ; overlay(plot_Qh3, plot_Qle3)

  plots(0) = plot_Qh1
  plots(1) = plot_Qh2
  ; plots(2) = plot_Qh3

else if LH_SH_scale .eq. "monthly" then
  time_tot        = year_sum*12
  res1@tmXBValues  = ispan(1,time_tot+1,12*2)
  x               = ispan(1,time_tot,1)
  res1@tmXBLabels  = ispan(year_s,year_e+1,2)
  ; res1@tiXAxisString = "Year"
  plots(0)          = gsn_csm_xy(wks, x, Qle_season(1,0,:)-Qle_season(0,0,:),res1)
  plot1             = gsn_csm_xy(wks, x, Qh_season(1,0,:)-Qh_season(0,0,:),res1)
  overlay(plots(0), plot1)
  plots(1)          = gsn_csm_xy(wks, x, Qle_season(1,1,:)-Qle_season(0,1,:),res1)
  plot2             = gsn_csm_xy(wks, x, Qh_season(1,1,:)-Qh_season(0,1,:),res1)
  overlay(plots(1), plot2)
end if
end if

; _____________________________ plot EVAPFBL line ______________________________

if EVAPFBL_plot .eq. "line" then

  res2@xyLineColors      = (/"red","orange","pink","green","blue","purple"/);,\
                            ; "red","orange","pink","green","blue","purple"/)
  res2@xyLineThicknesses = (/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/);,\
                            ;  1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
  res2@xyDashPatterns    = (/ 0,   0,   0,   0,   0,   0/)
                          ;(/   0,   0,   0,   0,   0,   0,\
                          ;      1,   1,   1,   1,   1,   1/)
  res2@xyLineOpacities   = (/1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
                           ;(/ 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,\
                            ;  1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)

  textres                   = True
  gres                      = True
  lineres                   = True
  lineres@lgLineColors      = (/"red","orange","pink","green","blue","purple"/)
  lineres@lgLineThicknesses = 1.    ; line thicknesses
  lineres@LineLengthPercent = 3.    ; expressed as %, 0->100, length of line
  gres@YPosPercent          = 90.  ;expressed as %, 0->100, sets position of top border of legend (Default = 95.)
  gres@XPosPercent          = 50.  ;expressed as %, 0->100, sets position of left border of legend (Default = 5.)
  lineres@lgDashIndexes     = (/0,0,0,0,0,0/)
  textres@lgLabels          = layer_depth

  ; res2@lgPerimOn              = False      ; No legend perimeter.
  ; res2@pmLegendParallelPosF   = 0.25
  ; res2@pmLegendOrthogonalPosF = -1.2
  ; res2@lgBoxMinorExtentF      = 0.2
  ; res2@pmLegendWidthF         = 0.15
  ; res2@pmLegendHeightF        = 0.15
  ; res2@lgItemOrder            = (/5,4,3,2,1,0/)   ; Reorder the legends
  ; res2@xyExplicitLegendLabels = layer_depth

  ; res2@pmLegendDisplayMode = "Always"
  res2@trYMinF             = -10.
  res2@trYMaxF             = 16.

  if EVAPFBL_scale .eq. "seasonal" then
    season_tot       = year_sum*4
    x                = ispan(1,season_tot,1)-1./5.
    res2@tmXBValues  = ispan(1,season_tot+1,4*2)
    res2@tmXBLabels  = ispan(year_s,year_e+1,2)
    res2@trXMinF     = 1.
    res2@trXMaxF     = season_tot+1

    ; plots(3)         = gsn_csm_xy(wks, x, EVAPFBL_season(0,1,:,:), res2)
    ; plots(4)         = gsn_csm_xy(wks, x, EVAPFBL_season(1,1,:,:), res2)
    ; plots(5)         = gsn_csm_xy(wks, x, EVAPFBL_season(2,1,:,:), res2)
    ; plots(6)         = gsn_csm_xy(wks, x, EVAPFBL_season(3,1,:,:), res2)

    ; FD - CTL
    ; EVAPFBL2  = new((/12,year_sum*4/),float, default_val)
    ; EVAPFBL2(0:5,:)  = EVAPFBL_season(0,1,:,:) ; CTL
    ; EVAPFBL2(6:11,:) = EVAPFBL_season(1,1,:,:) ; FD
    res2@tiYAxisString = "~F5~D"+"~F25~E (mm season~S~-1~N~)"
    plots(2)         = gsn_csm_xy(wks, x, EVAPFBL_season(1,1,:,:) - EVAPFBL_season(0,1,:,:), res2)
    plots(2)         = simple_legend(wks,plots(2),gres,lineres,textres)
    delete(res2@tiYAxisString)
    ; delete(EVAPFBL2)

    ; DR - CTL
    res2@pmLegendDisplayMode = "Never"
    ; EVAPFBL2  = new((/12,year_sum*4/),float, default_val)
    ; EVAPFBL2(0:5,:)  = EVAPFBL_season(0,1,:,:) ; CTL
    ; EVAPFBL2(6:11,:) = EVAPFBL_season(2,1,:,:) ; DR
    plots(3)         = gsn_csm_xy(wks, x, EVAPFBL_season(2,1,:,:) - EVAPFBL_season(0,1,:,:), res2)
    ; delete(EVAPFBL2)

    ; FD+DR - FD
    ; EVAPFBL2  = new((/12,year_sum*4/),float, default_val)
    ; EVAPFBL2(0:5,:)  = EVAPFBL_season(1,1,:,:) ; FD
    ; EVAPFBL2(6:11,:) = EVAPFBL_season(3,1,:,:) ; FD+DR
    ; plots(5)         = gsn_csm_xy(wks, x, EVAPFBL_season(3,1,:,:) - EVAPFBL_season(1,1,:,:), res2)
    ; delete(EVAPFBL2)

  else if EVAPFBL_scale .eq. "monthly" then
    month_tot        = year_sum*12
    x                = ispan(1,month_tot,1)
    res2@tmXBValues  = ispan(1,month_tot+1,12*2)
    res2@tmXBLabels  = ispan(year_s,year_e+1,2)

    ; plots(3)         = gsn_csm_xy(wks, x, EVAPFBL(0,1,:,:), res2)
    ; plots(4)         = gsn_csm_xy(wks, x, EVAPFBL(1,1,:,:), res2)
    ; plots(5)         = gsn_csm_xy(wks, x, EVAPFBL(2,1,:,:), res2)
    ; plots(6)         = gsn_csm_xy(wks, x, EVAPFBL(3,1,:,:), res2)

    ; FD - CTL
    res2@tiYAxisString = "~F25~E (mm mth~S~-1~N~)";"~F5~D"+
    plots(2)         = gsn_csm_xy(wks, x, EVAPFBL(1,1,:,:) - EVAPFBL(0,1,:,:), res2)
    delete(res2@tiYAxisString)

    ; DR - CTL
    res2@pmLegendDisplayMode = "Never"
    plots(3)         = gsn_csm_xy(wks, x, EVAPFBL(2,1,:,:) - EVAPFBL(0,1,:,:), res2)
    ; FD+DR - FD
    ; plots(5)         = gsn_csm_xy(wks, x, EVAPFBL(3,1,:,:) - EVAPFBL(1,1,:,:), res2)

  end if
  end if

else if EVAPFBL_plot .eq. "bar" then

  ; res@vpXF           = 0.03
  ; res@vpYF           = 0.7
  res2@gsnYRefLine    = 0.              ; reference line
  res2@gsnXYBarChart  = True            ; create bar chart
  res2@trYMinF        = 0.
  res2@trYMaxF        = 600.
  res2@pmLegendDisplayMode  = "Never"

  delete(x)
  year_tot         = year_sum*2
  res2@tmXBValues  = ispan(1,year_tot+1,2*2)
  x                = ispan(1,year_tot,1)
  res2@tmXBLabels  = ispan(year_s,year_e+1,2)
  res2@trXMinF     = 1.
  res2@trXMaxF     = year_tot+1


  ; All layers
  res2@gsnXYBarChartColors2 = (/"blue","forestgreen"/)
  res2@xyLineColors         = (/"white","forestgreen"/)
  res2@tiYAxisString        = "Et (mm yr~S~-1~N~)"

  res2@gsnXYBarChartFillOpacityF      = 0.3
  plot_all_layers1          = gsn_csm_xy(wks, x, EVAPFBL_bar(1,0,0,:),res2) ; tree, CTL&FD, all layers
  delete(res2@tiYAxisString)
  plot_all_layers2          = gsn_csm_xy(wks, x, EVAPFBL_bar(1,1,0,:),res2) ; tree, CTL&DR, all layers

  ; Deep soils
  res2@gsnXYBarChartFillOpacityF      = 1.
  plot_deep_soil1           = gsn_csm_xy(wks, x, EVAPFBL_bar(1,0,1,:),res2) ; tree, CTL&FD, all layers
  plot_deep_soil2           = gsn_csm_xy(wks, x, EVAPFBL_bar(1,1,1,:),res2) ; tree, CTL&DR, all layers

  overlay(plot_deep_soil1, plot_all_layers1)
  overlay(plot_deep_soil2, plot_all_layers2)

  plots(2) = plot_deep_soil1
  plots(3) = plot_deep_soil2

  textres                   = True
  gres                      = True
  lineres                   = True
  lineres@lgLineColors      = (/"blue","blue","forestgreen","forestgreen"/)
  lineres@lgLineThicknesses = 10.    ; line thicknesses
  lineres@lgLabelFontHeightF = 0.1
  lineres@LineLengthPercent = 2.    ; expressed as %, 0->100, length of line
  gres@YPosPercent          = 95.  ;expressed as %, 0->100, sets position of top border of legend (Default = 95.)
  gres@XPosPercent          = 30.  ;expressed as %, 0->100, sets position of left border of legend (Default = 5.)
  lineres@lgDashIndexes     = (/0,0,0,0/)
  textres@lgLabels          = (/"CTL Shollow","CTL Deep","FD Shollow","FD Deep"/)
  plots(2) = simple_legend(wks,plots(2),gres,lineres,textres)
  textres@lgLabels          = (/"CTL Shollow","CTL Deep","DR Shollow","DR Deep"/)
  plots(3) = simple_legend(wks,plots(3),gres,lineres,textres)
  delete(x)
else if EVAPFBL_plot .eq. "profile" then

  ; ___________________________ Plot EVAPFBL profile _____________________________
  ; plot EVAPFBL diff -- TREE
  res@lbTitleString  = "T ~C~(mm d~S~-1~N~)"
  res@cnLevels       = fspan(-0.1,0.1,21)
  res@cnFillPalette  = "MPL_BrBG" ;"ViBlGrWhYeOrRe"
  res@tiYAxisString  = "Depth (cm)"
  ; res@tiXAxisString  = "Year"

  ; FD - CTL
  res@lbLabelBarOn   = False
  tmp                = EVAPFBL_regrid(1,1,::-1,:) - EVAPFBL_regrid(0,1,::-1,:)
  plots(6)           = gsn_csm_contour(wks, tmp*1., res)
  delete(res@tiYAxisString)
  delete(tmp)

  ; DR - CTL
  tmp                = EVAPFBL_regrid(2,1,::-1,:) - EVAPFBL_regrid(0,1,::-1,:)
  plots(7)           = gsn_csm_contour(wks, tmp*1., res)
  delete(tmp)

  ; ; FD+DR - FD
  ; res@lbLabelBarOn   = True
  ; tmp                = EVAPFBL_regrid(3,1,::-1,:) - EVAPFBL_regrid(1,1,::-1,:)
  ; plots(8)           = gsn_csm_contour(wks, tmp*1., res)
  ; delete(tmp)

end if
end if
end if
; ______________________________ Plot SM profile _______________________________
; ; plot SM
; res@cnLevels       = fspan(0.0,0.4,21)
; res@cnFillPalette  = "GMT_drywet" ;"WhiteBlueGreenYellowRed"
; res@tiYAxisString  = "Depth (cm)"
; ; res@lbTitleString  = "(m3/m3)"
; plots(2)           = gsn_csm_contour(wks, SoilMoist_regrid(0,0,::-1,:)*1.,res)
; delete(res@tiYAxisString)
; plots(3)           = gsn_csm_contour(wks, SoilMoist_regrid(0,1,::-1,:)*1.,res)
; delete(res@cnLevels)
; ; plots(2)           = gsn_csm_contour(wks, SoilMoist_regrid(1,0,::-1,:)*1.,res)
; ; plots(3)           = gsn_csm_contour(wks, SoilMoist_regrid(1,1,::-1,:)*1.,res)

; plot SM diff -- TREE
res@lbTitleString  = "~F5~D"+"~F25~SM ~C~(m~S~3~N~ m~S~-3~N~)"
res@cnLevels       = fspan(-0.1,0.1,21)
res@cnFillPalette  = "MPL_BrBG" ;"ViBlGrWhYeOrRe"
; res@tiXAxisString  = "Year"

res@lbLabelBarOn   = False

; FD - CTL
res@tiYAxisString  = "Depth (cm)"
tmp                = SoilMoist_regrid(1,1,::-1,:) - SoilMoist_regrid(0,1,::-1,:)
plots(4)           = gsn_csm_contour(wks, tmp*1., res)
delete(res@tiYAxisString)
delete(tmp)

; DR - CTL
tmp                = SoilMoist_regrid(2,1,::-1,:) - SoilMoist_regrid(0,1,::-1,:)
plots(5)           = gsn_csm_contour(wks, tmp*1., res)
delete(tmp)

; ; FD+DR - FD
; res@lbLabelBarOn   = True
; tmp                = SoilMoist_regrid(3,1,::-1,:) - SoilMoist_regrid(1,1,::-1,:)
; plots(8)           = gsn_csm_contour(wks, tmp*1., res)
; delete(tmp)

do counter = 0,5
  txid1   = gsn_create_text(wks, order_name(counter), txres)
  annoid1 = gsn_add_annotation(plots(counter), txid1, amres)
end do

gsn_panel(wks,(/plots/),(/3,2/),pres)
exit
end
