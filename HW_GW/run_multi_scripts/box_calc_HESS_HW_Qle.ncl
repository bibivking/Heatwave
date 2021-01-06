;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./boxplot_mmy.ncl"
load "./boxplot_stats_mmy.ncl"
;*******************************************

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; ================================== Options ===================================
year_s       = 2000
year_e       = 2019
operation    = "calc"    ;"calc" ; "plot"
scale        = "SE-AU"   ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
heatwave     = "hw-90th" ; "hw-90th" ; "all-day"
case_name    = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr", \
                 "CTL-30x3+30yr-DeepRoot/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
var_name     = "Qle"

iveg_tree_s = 1
iveg_tree_e = 4
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

; _________ 4 time periods __________
          ;  all years,  2001-2009,  2017-2019,   2010-2011
          ;   2000-1-1,   2001-1-1,   2017-1-1,    2010-1-1
ts_s    = (/         0,        366,       6210,       3653/)
          ; 2019-12-31, 2009-12-31, 2019-12-31,  2011-12-31
ts_e    = (/ day_sum-1,       3652,       7304,       4382/)


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
filename  = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f         = addfile (filename,"r")
lat_sum   = dimsizes(f->Rainf(0,{A:B},0))
lon_sum   = dimsizes(f->Rainf(0,0,{C:D}))
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

total_box      = case_sum*2*4+3 ; case_sum x PFT all/tree x 4 time periods + space

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

; _____________________________ Calc Box Stats ________________________________
var_stats      = new((/total_box, 5/), float, -9999.)
var_stats      = -9999.

do case_num = 0,case_sum-1

  print("case_num = " +case_num)
  ; ______________________________ SET DATA _________________________________
  var_val    = new((/day_sum,lat_sum,lon_sum/),float,def_value)
  var_val    = 0.
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

    if (var_name .eq. "Qle") .or. (var_name .eq. "Qh") .or. (var_name .eq. "Fwsoil") then
      var_val(day_step_s:day_step_e,:,:)   = f1->$var_name$(:,{A:B},{C:D})
    else if var_name .eq. "deltaT" then
      var_val(day_step_s:day_step_e,:,:)   = f1->VegT(:,{A:B},{C:D}) - f1->Tair(:,{A:B},{C:D})
    else if var_name .eq. "SM64cm" then
      var_val(day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022  \
                                              + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
                                              + f1->SoilMoist(:,2,{A:B},{C:D})*0.154 \
                                              + f1->SoilMoist(:,3,{A:B},{C:D})*0.409)\
                                              / 0.643
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
    print("mask HW")
    var_val      = where(HW .eq. 1, var_val, def_value)
  end if

  ; _______________________________ MASK PFT __________________________________
  print("mask PFT")
  ; trees
  var_val_tree      = where(PFT_mask_multi .eq. 1, var_val, def_value)

  ; non-trees
  var_val      = where(PFT_mask_multi .eq. 2, var_val, def_value)

  loc1 = case_num
  loc2 = case_num + case_sum
  groups = case_sum*2+1

  do count = 0, 3
    var_stats(count*groups+loc1, :)    = boxplot_stats_mmy(var_val(ts_s(count):ts_e(count),:,:),def_value)
    var_stats(count*groups+loc2, :)    = boxplot_stats_mmy(var_val_tree(ts_s(count):ts_e(count),:,:),def_value)
  end do ; count

  delete(var_val)
  delete(var_val_tree)
end do ; case_num

if heatwave .eq. "hw-90th" then
  delete(HW)
end if
delete(PFT_mask_multi)

; _____________________________ Write box stats ________________________________
print("write out")
option        = True
option@row    = False
; option@tspace = 0

var_stats   = where( ismissing(var_stats), -9999., var_stats)
print(var_stats)
var_stats@_FillValue = -9999.
option@fout = var_name+"_box_stats.txt"
format      = case_sum+"f20.4"
write_matrix (var_stats, format, option)
delete(option@fout)

end
