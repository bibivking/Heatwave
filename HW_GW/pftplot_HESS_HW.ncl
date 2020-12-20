;*****************************************
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

undef("pdf_1d_mmy")
function pdf_1d_mmy(var[*][*][*]:numeric,bins:numeric,def_value:numeric)
;
;    var[*][*][*] : input data
;    bins         : total bin numbers
;    def_value    : default values
;
begin

  var@_FillValue = def_value

  ; ___________ nd to 1d ___________
  tmp            = ndtooned(var)
  maskout_def    = ind( .not. ismissing(tmp))
  print(dimsizes(maskout_def))

  var_masked     = tmp(maskout_def)
  var_masked@_FillValue = def_value

  ; ___________ pdf ___________
  var_pdf = pdfx(var_masked, bins, False)

  return(var_pdf)
end


; ================================== Options ===================================
year_s       = 2000
year_e       = 2019

scale        = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
heatwave     = "hw-90th" ; "hw-90th" ; "all-day"
case_name    = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr"/);,"PUMPx10","PUMPx100"/)
message      = "CTL_vs_FREE_DRAIN"
pic_name     = "HW_pdfplot_all_"+heatwave+"_"+scale+"_"+message+"_"+year_s+"-"+year_e

iveg_s = 1
iveg_e = 4

; ================================= Predefine ==================================
case_sum  = dimsizes(case_name)
year_sum = year_e-year_s+1

day_sum  = 0
do year = year_s,year_e
  if ( mod(year,4) .eq. 0) then
     day_sum = day_sum + 366
  else
     day_sum = day_sum + 365
  end if
end do


; =================================== Path =====================================
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"+case_name(case_num)
end do

; ================================ Spitial Scale ===============================
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

; ================================== Start =====================================
filename = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f        = addfile (filename,"r")
lat_sum  = dimsizes(f->Rainf(0,{A:B},0))
lon_sum  = dimsizes(f->Rainf(0,0,{C:D}))
def_value = f->Rainf@_FillValue

; PFT mask
PFT_mask       = where((f->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
PFT_mask_multi = new((/day_sum,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, day_sum-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do
delete(PFT_mask)
delete(f)


; ______________________________ READ HW _________________________________
if heatwave .eq. "hw-90th" then
  filename_hw = "/g/data/w35/mm3972/scripts/ehfheatwaves/HW_Event_Indicator_1970-2019.nc" ; summer HW only
  fhw         = addfile(filename_hw,"r")

  HW = new((/day_sum,lat_sum,lon_sum/),float)
  day_start = 0
  do year = 1970, year_s - 1
    if ( mod(year,4) .eq. 0) then
       day_start = day_start + 366
    else
       day_start = day_start + 365
    end if
  end do
  HW = fhw->event(day_start:day_start+day_sum-1,:,:)*1.
  delete(filename_hw)
  delete(fhw)
end if

; _____________________________ Calc Bin Stats ________________________________
bins     = 40
nVar     = 16
nBin     = bins

Qle_x      = new((/nVar, nBin/),integer)
Qh_x       = new((/nVar, nBin/),integer)
deltaT_x   = new((/nVar, nBin/),integer)
Fwsoil_x   = new((/nVar, nBin/),integer)
SM64cm_x   = new((/nVar, nBin/),integer)
WatTable_x = new((/nVar, nBin/),integer)

Qle_y      = new((/nVar, nBin/),float,def_value)
Qh_y       = new((/nVar, nBin/),float,def_value)
deltaT_y   = new((/nVar, nBin/),float,def_value)
Fwsoil_y   = new((/nVar, nBin/),float,def_value)
SM64cm_y   = new((/nVar, nBin/),float,def_value)
WatTable_y = new((/nVar, nBin/),float,def_value)

          ;  all years,  2001-2009,  2017-2019,   2010-2011
          ;   2000-1-1,   2001-1-1,   2017-1-1,    2010-1-1
ts_s    = (/         0,        366,       6210,       3653/)
          ; 2019-12-31, 2009-12-31, 2019-12-31,  2011-12-31
ts_e    = (/ day_sum-1,       3652,       7304,       4382/)


; _____________________________ Start Main Loop ________________________________
do case_num = 0,case_sum-1

  print("case_num = " +case_num)

  ; ______________________________ SET DATA _________________________________
  Qle      = new((/day_sum,lat_sum,lon_sum/),float,def_value)
  Qh       = new((/day_sum,lat_sum,lon_sum/),float,def_value)
  deltaT   = new((/day_sum,lat_sum,lon_sum/),float,def_value)
  Fwsoil   = new((/day_sum,lat_sum,lon_sum/),float,def_value)
  SM64cm   = new((/day_sum,lat_sum,lon_sum/),float,def_value)
  WatTable = new((/day_sum,lat_sum,lon_sum/),float,def_value)

  Qle      = 0.
  Qh       = 0.
  deltaT   = 0.
  Fwsoil   = 0.
  SM64cm   = 0.
  WatTable = 0.

  day_step_s = 0

  ; ______________________________ READ DATA _________________________________
  do year = year_s,year_e
    print("year = " +year)

    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    ; calc last day of year
    yr_tot_day = sum(dom)
    day_step_e = day_step_s + yr_tot_day - 1

    filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
    f1 = addfile (filename1,"r")

    Qle(day_step_s:day_step_e,:,:)      = f1->Qle(:,{A:B},{C:D})
    Qh(day_step_s:day_step_e,:,:)       = f1->Qh(:,{A:B},{C:D})
    Fwsoil(day_step_s:day_step_e,:,:)   = f1->Fwsoil(:,{A:B},{C:D})
    deltaT(day_step_s:day_step_e,:,:)   = f1->VegT(:,{A:B},{C:D}) - f1->Tair(:,{A:B},{C:D})
    WatTable(day_step_s:day_step_e,:,:) = f1->WatTable(:,{A:B},{C:D})
    SM64cm(day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022  \
                                                    + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
                                                    + f1->SoilMoist(:,2,{A:B},{C:D})*0.154 \
                                                    + f1->SoilMoist(:,3,{A:B},{C:D})*0.409)\
                                                    / 0.643
    ; calc fist day of next year
    day_step_s = day_step_e + 1

    delete(filename1)
    delete(f1)
  end do ; year

  ; _____________________________ MASK HW DAY ________________________________
  if heatwave .eq. "hw-90th" then
    Qle      = where(HW .eq. 1, Qle, def_value)
    Qh       = where(HW .eq. 1, Qh, def_value)
    Fwsoil   = where(HW .eq. 1, Fwsoil, def_value)
    deltaT   = where(HW .eq. 1, deltaT, def_value)
    WatTable = where(HW .eq. 1, WatTable, def_value)
    SM64cm   = where(HW .eq. 1, SM64cm, def_value)
  end if

  ; _______________________________ MASK PFT __________________________________
  Qle_tree      = where(PFT_mask_multi .eq. 1, Qle, def_value)
  Qh_tree       = where(PFT_mask_multi .eq. 1, Qh, def_value)
  Fwsoil_tree   = where(PFT_mask_multi .eq. 1, Fwsoil, def_value)
  deltaT_tree   = where(PFT_mask_multi .eq. 1, deltaT, def_value)
  WatTable_tree = where(PFT_mask_multi .eq. 1, WatTable, def_value)
  SM64cm_tree   = where(PFT_mask_multi .eq. 1, SM64cm, def_value)

  if case_num .eq. 0 then
    loc1 = 0
    loc2 = 2
  else if case_num .eq. 1 then
    loc1 = 1
    loc2 = 3
  end if
  end if

  do count = 0, 3 ; four periods

    tmp = pdf_1d_mmy(Qle(ts_s(count):ts_e(count),:,:),bins, def_value)
    print(tmp@bin_center)
    Qle_x(count*4+loc1,:) = tmp@bin_center
    Qle_y(count*4+loc1,:) = tmp
    delete(tmp)
    tmp = pdf_1d_mmy(Qle_tree(ts_s(count):ts_e(count),:,:),bins, def_value)
    Qle_x(count*4+loc2,:) = tmp@bin_center
    Qle_y(count*4+loc2,:) = tmp
    delete(tmp)

    tmp = pdf_1d_mmy(Qh(ts_s(count):ts_e(count),:,:),bins, def_value)
    Qh_x(count*4+loc1,:) = tmp@bin_center
    Qh_y(count*4+loc1,:) = tmp
    delete(tmp)
    tmp = pdf_1d_mmy(Qh_tree(ts_s(count):ts_e(count),:,:),bins, def_value)
    Qh_x(count*4+loc2,:) = tmp@bin_center
    Qh_y(count*4+loc2,:) = tmp
    delete(tmp)

    tmp = pdf_1d_mmy(Fwsoil(ts_s(count):ts_e(count),:,:),bins, def_value)
    Fwsoil_x(count*4+loc1,:) = tmp@bin_center
    Fwsoil_y(count*4+loc1,:) = tmp
    delete(tmp)
    tmp = pdf_1d_mmy(Fwsoil_tree(ts_s(count):ts_e(count),:,:),bins, def_value)
    Fwsoil_x(count*4+loc2,:) = tmp@bin_center
    Fwsoil_y(count*4+loc2,:) = tmp
    delete(tmp)

    tmp = pdf_1d_mmy(deltaT(ts_s(count):ts_e(count),:,:),bins, def_value)
    deltaT_x(count*4+loc1,:) = tmp@bin_center
    deltaT_y(count*4+loc1,:) = tmp
    delete(tmp)
    tmp = pdf_1d_mmy(deltaT_tree(ts_s(count):ts_e(count),:,:),bins, def_value)
    deltaT_x(count*4+loc2,:) = tmp@bin_center
    deltaT_y(count*4+loc2,:) = tmp
    delete(tmp)

    tmp = pdf_1d_mmy(WatTable(ts_s(count):ts_e(count),:,:),bins, def_value)
    WatTable_x(count*4+loc1,:) = tmp@bin_center
    WatTable_y(count*4+loc1,:) = tmp
    delete(tmp)
    tmp = pdf_1d_mmy(WatTable_tree(ts_s(count):ts_e(count),:,:),bins, def_value)
    WatTable_x(count*4+loc2,:) = tmp@bin_center
    WatTable_y(count*4+loc2,:) = tmp
    delete(tmp)

    tmp = pdf_1d_mmy(SM64cm(ts_s(count):ts_e(count),:,:),bins, def_value)
    SM64cm_x(count*4+loc1,:) = tmp@bin_center
    SM64cm_y(count*4+loc1,:) = tmp
    delete(tmp)
    tmp = pdf_1d_mmy(SM64cm_tree(ts_s(count):ts_e(count),:,:),bins, def_value)
    SM64cm_x(count*4+loc2,:) = tmp@bin_center
    SM64cm_y(count*4+loc2,:) = tmp
    delete(tmp)

  end do ; count

  delete(Qle)
  delete(Qh)
  delete(deltaT)
  delete(Fwsoil)
  delete(SM64cm)
  delete(WatTable)

  delete(Qle_tree)
  delete(Qh_tree)
  delete(Fwsoil_tree)
  delete(deltaT_tree)
  delete(WatTable_tree)
  delete(SM64cm_tree)
end do ; case_num

if heatwave .eq. "hw-90th" then
  delete(HW)
end if
delete(PFT_mask_multi)


;_____________________________ PLOT ___________________________________
wks = gsn_open_wks("pdf",pic_name)
gsn_define_colormap(wks,"rainbow") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")

res                          = True
res@gsnDraw                  = True  ; Don't draw plots
res@gsnFrame                 = True  ; ����ҳ
res@gsnXYBarChart            = True              ; Create bar plot
res@gsnXYBarChartOutlineOnly = True

res@pmLegendDisplayMode    = "Always"            ; turn on legend
res@pmLegendSide           = "Top"               ; Change location of
res@pmLegendParallelPosF   = .25                 ; move units right
res@pmLegendOrthogonalPosF = -0.4                ; move units down
res@pmLegendWidthF         = 0.125               ; Change width and
res@pmLegendHeightF        = 0.15                ; height of legend.
res@lgPerimOn              = True                ; turn off/on box around
res@lgLabelFontHeightF     = .015                ; label font height
res@tiYAxisString          = "PDF (%)"
res@tiMainString           = "PDF: Three Distributions: 40 Bins"

res@xyLineThicknesses      = (/2.0,2.0,2.0,2.0,\
                               2.0,2.0,2.0,2.0,\
                               2.0,2.0,2.0,2.0,\
                               2.0,2.0,2.0,2.0 /)
res@xyLineColors           = (/"blue","red","green","orange",\
                               "blue","red","green","orange",\
                               "blue","red","green","orange",\
                               "blue","red","green","orange" /)
res@xyDashPatterns         = (/0,0,0,0,\
                               1,1,1,1,\
                               2,2,2,2,\
                               3,3,3,3 /)

res@xyExplicitLegendLabels = (/"CTL_ALL_All","FD_ALL_All","CTL_TREE_All","FD_TREE_All",\
                              "CTL_ALL_01-09","FD_ALL_01-09","CTL_TREE_01-09","FD_TREE_01-09",\
                              "CTL_ALL_17-19","FD_ALL_17-19","CTL_TREE_17-19","FD_TREE_17-19",\
                              "CTL_ALL_10-11","FD_ALL_10-11","CTL_TREE_10-11","FD_TREE_10-11"/)  ; create explicit labels

plot = gsn_csm_xy (wks, Qle_x, Qle_y, res)
plot = gsn_csm_xy (wks, Qh_x, Qh_y, res)
plot = gsn_csm_xy (wks, deltaT_x, deltaT_y, res)
plot = gsn_csm_xy (wks, Fwsoil_x, Fwsoil_y, res)
plot = gsn_csm_xy (wks, SM64cm_x, SM64cm_y, res)
plot = gsn_csm_xy (wks, WatTable_x, WatTable_y, res)

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
exit
end
