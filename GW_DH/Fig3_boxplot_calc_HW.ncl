
;MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
; Author  : Mengyuan Mu
; Version : 2.0 (23.06.2021)"
; Email   : mu.mengyuan815@gmail.com
;WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./Fig3_boxplot_stats_mmy.ncl"
;*******************************************

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; ================================== Options ===================================
year_s       = 2000
year_e       = 2001;19
operation    = "calc"    ;"calc" ; "plot"
scale        = "SE-AU"   ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
heatwave     = "hw-90th" ; "hw-90th" ; "all-day"
case_name    = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr"/)
case_short_name = (/"GW","FD"/)
var_name     = "EF" ;"TVeg" ; "GPP"
select_region= False
median_diff = True
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
         ; 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,\
         ; 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019
ts_s    = (/335, 700, 1065, 1430, 1796, 2161, 2526, 2891, 3257, 3622,\
           3987,4352, 4718, 5083, 5448, 5813, 6179, 6544, 6909/) ; ;

ts_e    = (/425, 790, 1155, 1521, 1886, 2251, 2616, 2982, 3347, 3712, \
           4077, 4443, 4808,5173, 5538, 5904, 6269, 6634, 6999/) ; ;


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
def_iveg  = f->iveg@_FillValue

; PFT mask
if select_region .eq. True then
  LAI      = dim_avg_n_Wrap(f->LAI(:,{A:B},{C:D}),0)
  PFT_mask = where(LAI .gt. 1.5, f->iveg({A:B},{C:D}), def_iveg)
  PFT_mask = where((f->latitude({A:B},{C:D}) .le. -35.) .and. (f->latitude({A:B},{C:D}) .ge. -40.)   \
            .and.  (f->longitude({A:B},{C:D}) .le. 149.) .and. (f->longitude({A:B},{C:D}) .ge. 145.),\
             PFT_mask, def_iveg)
  PFT_mask = where((PFT_mask .ge. iveg_tree_s) .and. (PFT_mask .le. iveg_tree_e), 1, PFT_mask)
  PFT_mask = where((PFT_mask .ge. iveg_nontree_s) .and. (PFT_mask .le. iveg_nontree_e), 2, PFT_mask)
else
  PFT_mask = where((f->iveg({A:B},{C:D}) .ge. iveg_tree_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_tree_e), 1, f->iveg({A:B},{C:D}))
  PFT_mask = where((PFT_mask .ge. iveg_nontree_s) .and. (PFT_mask .le. iveg_nontree_e), 2, PFT_mask)
end if

PFT_mask_multi = new((/day_sum,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, day_sum-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do
; delete(PFT_mask)
delete(f)

total_box      = case_sum*2*dimsizes(ts_s)+dimsizes(ts_s)-1 ;case_sum*2*4+3 ; case_sum x PFT all/tree x 4 time periods + space

; ______________________________ READ HW _________________________________
if heatwave .eq. "hw-90th" then
  filename_hw = "/g/data/w35/mm3972/scripts/ehfheatwaves/nc_file/SE_AUS/HW_Event_Indicator_1970-2019.nc" ; summer HW only
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

; ====================== Checking mask ===========================
wks = gsn_open_wks("pdf","Check_boxplot")
gsn_define_colormap(wks,"rainbow")

res                         = True
res@gsnDraw                 = True  ; Don't draw plots
res@gsnFrame                = True  ; ����ҳ
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

plot = gsn_csm_contour(wks,PFT_mask,res)

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
    end if
    if (var_name .eq. "Evap") .or. (var_name .eq. "TVeg") .or. (var_name .eq. "ESoil") then
      ET_scale = 3600*24. ; /s ->/day
      var_val(day_step_s:day_step_e,:,:)   = f1->$var_name$(:,{A:B},{C:D})*ET_scale
    end if
    if (var_name .eq. "GPP") .or. (var_name .eq. "NPP") then
      s2d       = 3600*24. ; s-1 to d-1
      GPP_scale = 0.000001*12*s2d ; umol s-1 to g d-1
      var_val(day_step_s:day_step_e,:,:)   = f1->$var_name$(:,{A:B},{C:D})*GPP_scale
    end if
    if var_name .eq. "deltaT" then
      var_val(day_step_s:day_step_e,:,:)   = f1->VegT(:,{A:B},{C:D}) - f1->Tair(:,{A:B},{C:D})
    end if
    if var_name .eq. "SM64cm" then
      var_val(day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022  \
                                              + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
                                              + f1->SoilMoist(:,2,{A:B},{C:D})*0.154 \
                                              + f1->SoilMoist(:,3,{A:B},{C:D})*0.409)\
                                              / 0.643
    end if
    if var_name .eq. "EF" then
      Qle_tmp = f1->Qle(:,{A:B},{C:D})
      Qh_tmp  = f1->Qh(:,{A:B},{C:D})
      var_val(day_step_s:day_step_e,:,:)   = where( (Qle_tmp +Qh_tmp) .ne. 0, \
                                                     Qle_tmp/(Qh_tmp+ Qle_tmp), \
                                                     def_value)
      delete(Qle_tmp)
      delete(Qh_tmp)
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
  var_val_tree = where(PFT_mask_multi .eq. 1, var_val, def_value)

  ; ; non-trees
  ; var_val      = where(PFT_mask_multi .eq. 2, var_val, def_value)

  ; loc1 = case_num
  ; loc2 = case_num + case_sum
  ; groups = case_sum*2+1

  ; do count = 0, dimsizes(ts_s)-1
  ;   print(count)
  ;   var_stats(count*groups+loc1, :)    = boxplot_stats_mmy(var_val(ts_s(count):ts_e(count),:,:),def_value)
  ;   var_stats(count*groups+loc2, :)    = boxplot_stats_mmy(var_val_tree(ts_s(count):ts_e(count),:,:),def_value)
  ; end do ; count
  delete(var_val)

  ; _____________________________ Write raw data for python script ________________________________
  print("write raw data for Python script")
 
  date_3D = var_val_tree
  do day  = 1, day_sum
    date_3D(day-1,:,:) = day
  end do
  
  var_tmp        = ndtooned(var_val_tree)
  delete(var_val_tree)
  date_tmp       = ndtooned(date_3D)
  print(date_tmp)
  var_tmp@_FillValue = def_value
  maskout_def    = ind( .not. ismissing(var_tmp))

  var_1D         = var_tmp(maskout_def)
  date_1D        = date_tmp(maskout_def)

  print(dimsizes(var_1D))
  print(var_1D(0))
  print(dimsizes(date_1D)) 
  print(date_1D(0))
  ; output list 

  output_filename = "./txt/"+var_name+"_"+case_short_name(case_num)+"_rawdata_4_Python.txt"

  system("rm "+ output_filename)

  alist = [/date_1D, var_1D/]

  write_table(output_filename, "w", alist, "%6i, %10.4f")

end do ; case_num

if heatwave .eq. "hw-90th" then
  delete(HW)
end if
delete(PFT_mask_multi)

; ; _____________________________ Write box stats ________________________________
; print("write out")
; option        = True
; option@row    = False
; ; option@tspace = 0

; var_stats   = where( ismissing(var_stats), -9999., var_stats)
; print(var_stats)
; var_stats@_FillValue = -9999.
; option@fout = "./txt/"+var_name+"_CTL_FD_yearly_box_stats.txt"
; format      = "5f20.4"
; write_matrix (var_stats, format, option)
; delete(option@fout)

exit
end
