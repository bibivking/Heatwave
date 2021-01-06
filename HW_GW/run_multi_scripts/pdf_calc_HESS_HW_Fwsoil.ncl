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
function pdf_1d_mmy(var[*][*][*]:numeric,bins:numeric,def_value:numeric,bin_min:numeric,bin_max:numeric)
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
  opt = True
  opt@bin_min  = bin_min
  opt@bin_max  = bin_max

  var_pdf = pdfx(var_masked, bins, opt)

  return(var_pdf)
end

; ================================== Options ===================================
year_s       = 2000
year_e       = 2019

scale        = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
heatwave     = "hw-90th" ; "hw-90th" ; "all-day"
case_name    = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr", \
                 "CTL-30x3+30yr-DeepRoot/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
message      = "CTL_vs_FD_vs_DR"
var_name     = "Fwsoil"
pic_name     = "HW_check_pdfplot_"+var_name+"_"+heatwave+"_"+scale+"_"+message

iveg_tree_s  = 1
iveg_tree_e  = 4
iveg_nontree_s = 5
iveg_nontree_e = 17

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
PFT_mask = where((f->iveg({A:B},{C:D}) .ge. iveg_tree_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_tree_e), 1, f->iveg({A:B},{C:D}))
PFT_mask = where((PFT_mask .ge. iveg_nontree_s) .and. (PFT_mask .le. iveg_nontree_e), 2, PFT_mask)
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
bins     = 100
nVar     = 32

var_x      = new(bins,double)
var_y      = new((/nVar, bins/),double,def_value)

          ;  all years,  2001-2009,  2010-2011,   2017-2019
          ;   2000-1-1,   2001-1-1,   2010-1-1,    2017-1-1
ts_s    = (/         0,        366,       3653,       6210/)
          ; 2019-12-31, 2009-12-31, 2011-12-31,  2019-12-31
ts_e    = (/ day_sum-1,       3652,       4382,       7304/)

; _____________________________ Start Main Loop ________________________________
do case_num = 0,case_sum-1

  print("case_num = " +case_num)

  ; ______________________________ SET DATA _________________________________
  var      = new((/day_sum,lat_sum,lon_sum/),float,def_value)
  var      = 0.
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

    if var_name .eq. "Qle" then
      bin_min  =   0.
      bin_max  = 300.
      var(day_step_s:day_step_e,:,:)   = f1->$var_name$(:,{A:B},{C:D})
    else if var_name .eq. "Qh" then
      bin_min  = -50.
      bin_max  = 350.
      var(day_step_s:day_step_e,:,:)   = f1->$var_name$(:,{A:B},{C:D})
    else if var_name .eq. "Fwsoil" then
      bin_min  = 0.
      bin_max  = 1.
      var(day_step_s:day_step_e,:,:)   = f1->$var_name$(:,{A:B},{C:D})
    else if var_name .eq. "deltaT" then
      bin_min  = -5.
      bin_max  = 10.
      var(day_step_s:day_step_e,:,:)   = f1->VegT(:,{A:B},{C:D}) - f1->Tair(:,{A:B},{C:D})
    else if var_name .eq. "SM64cm" then
      bin_min  = 0.05
      bin_max  = 0.45
      var(day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022  \
                                          + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
                                          + f1->SoilMoist(:,2,{A:B},{C:D})*0.154 \
                                          + f1->SoilMoist(:,3,{A:B},{C:D})*0.409)\
                                          / 0.643
    end if
    end if
    end if
    end if
    end if

    ; calc fist day of next year
    day_step_s = day_step_e + 1

    delete(filename1)
    delete(f1)
  end do ; year

  ; _____________________________ MASK HW DAY ________________________________
  if heatwave .eq. "hw-90th" then
    print("mask hw-90th")
    var      = where(HW .eq. 1, var, def_value)
  end if

  ; _______________________________ MASK PFT __________________________________
  print("mask PFT")
  ; tree
  var_tree = where(PFT_mask_multi .eq. 1, var, def_value)
  ; non-trees
  var      = where(PFT_mask_multi .eq. 2, var, def_value)

  loc1 = case_num
  loc2 = case_num + case_sum

  groups = case_sum*2
  do count = 0, 3 ; four periods
    tmp = pdf_1d_mmy(var(ts_s(count):ts_e(count),:,:),bins,def_value,bin_min,bin_max)
    print(typeof(tmp@bin_center))

    if count .eq. 0 .and. case_num .eq. 0 then
      var_x                 = tmp@bin_center
    end if

    var_y(count*groups+loc1,:) = tmp
    delete(tmp)

    tmp = pdf_1d_mmy(var_tree(ts_s(count):ts_e(count),:,:),bins, def_value,bin_min,bin_max)
    var_y(count*groups+loc2,:) = tmp
    delete(tmp)
  end do ; count

  delete(var)
  delete(var_tree)
end do ; case_num

if heatwave .eq. "hw-90th" then
  delete(HW)
end if
delete(PFT_mask_multi)

; ______________________________ Write txt files _______________________________
print("write out")
option        = True
option@row    = False
; option@tspace = 0

var_y   = where( ismissing(var_y), -99999., var_y)
var_y@_FillValue = -99999.
var_all         = new((/nVar+1,bins/),double,def_value)
var_all(0,:)    = var_x
var_all(1:,:)   = var_y
option@fout = var_name + "_pdf.txt"
format      = bins+"f20.4"
write_matrix (var_all, format, option)
delete(option@fout)

;_____________________________ PLOT ___________________________________
wks = gsn_open_wks("pdf",pic_name)
gsn_define_colormap(wks,"rainbow") ;"BlueWhiteOrangeRed") ;"BlueYellowRed")

res                          = True
res@gsnDraw                  = True  ; Don't draw plots
res@gsnFrame                 = True  ; ����ҳ
; res@gsnXYBarChart            = True              ; Create bar plot
; res@gsnXYBarChartOutlineOnly = True

res@pmLegendDisplayMode    = "Always"            ; turn on legend
res@pmLegendSide           = "Right";"Top"               ; Change location of
res@pmLegendParallelPosF   = .25                 ; move units right
res@pmLegendOrthogonalPosF = -0.4                ; move units down
res@pmLegendWidthF         = 0.125               ; Change width and
res@pmLegendHeightF        = 0.15                ; height of legend.
res@lgPerimOn              = False                ; turn off/on box around
res@lgLabelFontHeightF     = .015                ; label font height
res@tiYAxisString          = "PDF (%)"

res@tmXBBorderOn      = False
res@tmXBOn            = False
res@tmXUseBottom      = False
res@tmXTLabelsOn      = True
; res@tmYLMinorPerMajor = 1
res@tmYRBorderOn      = False
res@tmYROn            = False
res@tmYUseLeft        = False
res@trYReverse        = True

res@xyLineThicknesses      = (/1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5/)
res@xyLineColors           = (/"blue","red","green","orange", "blue","red","green","orange"/)
res@xyDashPatterns         = (/1,1,1,1, 0,0,0,0/)

res@xyExplicitLegendLabels = (/"CTL_NonTREE","FD_NonTREE","DR_NonTREE","FD+DR_NonTREE", \
                               "CTL_TREE","FD_TREE","DR_TREE","FD+DR_TREE" /)  ; create explicit labels

res@tiMainString           = "2000-2019"
plot = gsn_csm_xy (wks, var_x, var_y(0:7,:), res)

res@tiMainString           = "2001-2009"
plot = gsn_csm_xy (wks, var_x, Qle_y(8:15,:), res)

res@tiMainString           = "2010-2011"
plot = gsn_csm_xy (wks, var_x, var_y(16:23,:), res)

res@tiMainString           = "2017-2019"
plot = gsn_csm_xy (wks, var_x, var_y(24:31,:), res)

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
