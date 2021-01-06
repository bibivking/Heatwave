;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./boxplot_mmy.ncl"
load "./boxplot_stats_mmy.ncl"
;*******************************************



; ================================= MMY =================================
; Aim : calculating statistic values for box plot 
;       and drawing boxplot              
; 
; =======================================================================

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; ================================== Options ===================================
year_s       = 2000
year_e       = 2019
operation    = "calc" ;"calc" ; "plot"
scale        = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
heatwave     = "hw-90th" ; "hw-90th" ; "all-day"
case_name    = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr", \
                 "CTL-30x3+30yr-DeepRoot/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
message      = "CTL_vs_FD_vs_DR"
pic_name     = "HW_boxplot_"+heatwave+"_"+scale+"_"+message+"_"+year_s+"-"+year_e

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

total_box      = case_sum*2*4+3 ; case_sum x PFT all/tree x 4 time periods + space

if operation .eq. "calc" then
  ; ============================ Calcute box stats =============================

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
  Qle_stats      = new((/total_box, 5/), float, -9999.)
  ; Qh_stats       = new((/total_box, 5/), float, -9999.)
  ; Fwsoil_stats   = new((/total_box, 5/), float, -9999.)
  ; deltaT_stats   = new((/total_box, 5/), float, -9999.)
  ; SM64cm_stats   = new((/total_box, 5/), float, -9999.)

  Qle_stats      = -9999.
  ; Qh_stats       = -9999.
  ; Fwsoil_stats   = -9999.
  ; deltaT_stats   = -9999.
  ; SM64cm_stats   = -9999.

  ; _____________________________ Start Main Loop ______________________________
  do case_num = 0,case_sum-1

    print("case_num = " +case_num)

    ; ______________________________ SET DATA _________________________________
    Qle      = new((/day_sum,lat_sum,lon_sum/),float,def_value)
    ; Qh       = new((/day_sum,lat_sum,lon_sum/),float,def_value)
    ; deltaT   = new((/day_sum,lat_sum,lon_sum/),float,def_value)
    ; Fwsoil   = new((/day_sum,lat_sum,lon_sum/),float,def_value)
    ; SM64cm   = new((/day_sum,lat_sum,lon_sum/),float,def_value)

    Qle      = 0.
    ; Qh       = 0.
    ; deltaT   = 0.
    ; Fwsoil   = 0.
    ; SM64cm   = 0.

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
      ; Qh(day_step_s:day_step_e,:,:)       = f1->Qh(:,{A:B},{C:D})
      ; Fwsoil(day_step_s:day_step_e,:,:)   = f1->Fwsoil(:,{A:B},{C:D})
      ; deltaT(day_step_s:day_step_e,:,:)   = f1->VegT(:,{A:B},{C:D}) - f1->Tair(:,{A:B},{C:D})
      ; SM64cm(day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022  \
      ;                                         + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
      ;                                         + f1->SoilMoist(:,2,{A:B},{C:D})*0.154 \
      ;                                         + f1->SoilMoist(:,3,{A:B},{C:D})*0.409)\
      ;                                         / 0.643
      ; calc fist day of next year
      day_step_s = day_step_e + 1

      delete(filename1)
      delete(f1)
    end do ; year

    ; _____________________________ MASK HW DAY ________________________________
    if heatwave .eq. "hw-90th" then
      print("mask HW")
      Qle      = where(HW .eq. 1, Qle, def_value)
      ; Qh       = where(HW .eq. 1, Qh, def_value)
      ; Fwsoil   = where(HW .eq. 1, Fwsoil, def_value)
      ; deltaT   = where(HW .eq. 1, deltaT, def_value)
      ; SM64cm   = where(HW .eq. 1, SM64cm, def_value)
    end if

    ; _______________________________ MASK PFT __________________________________
    print("mask PFT")
    ; trees
    Qle_tree      = where(PFT_mask_multi .eq. 1, Qle, def_value)
    ; Qh_tree       = where(PFT_mask_multi .eq. 1, Qh, def_value)
    ; Fwsoil_tree   = where(PFT_mask_multi .eq. 1, Fwsoil, def_value)
    ; deltaT_tree   = where(PFT_mask_multi .eq. 1, deltaT, def_value)
    ; SM64cm_tree   = where(PFT_mask_multi .eq. 1, SM64cm, def_value)

    ; non-trees
    Qle      = where(PFT_mask_multi .eq. 2, Qle, def_value)
    ; Qh       = where(PFT_mask_multi .eq. 2, Qh, def_value)
    ; Fwsoil   = where(PFT_mask_multi .eq. 2, Fwsoil, def_value)
    ; deltaT   = where(PFT_mask_multi .eq. 2, deltaT, def_value)
    ; SM64cm   = where(PFT_mask_multi .eq. 2, SM64cm, def_value)

    loc1 = case_num
    loc2 = case_num + case_sum
    groups = case_sum*2+1

    do count = 0, 3
      Qle_stats(count*groups+loc1, :)    = boxplot_stats_mmy(Qle(ts_s(count):ts_e(count),:,:),def_value)
      Qle_stats(count*groups+loc2, :)    = boxplot_stats_mmy(Qle_tree(ts_s(count):ts_e(count),:,:),def_value)

      ; Qh_stats(count*groups+loc1, :)     = boxplot_stats_mmy(Qh(ts_s(count):ts_e(count),:,:),def_value)
      ; Qh_stats(count*groups+loc2, :)     = boxplot_stats_mmy(Qh_tree(ts_s(count):ts_e(count),:,:),def_value)

      ; Fwsoil_stats(count*groups+loc1, :) = boxplot_stats_mmy(Fwsoil(ts_s(count):ts_e(count),:,:),def_value)
      ; Fwsoil_stats(count*groups+loc2, :) = boxplot_stats_mmy(Fwsoil_tree(ts_s(count):ts_e(count),:,:),def_value)
      ;
      ; deltaT_stats(count*groups+loc1, :) = boxplot_stats_mmy(deltaT(ts_s(count):ts_e(count),:,:),def_value)
      ; deltaT_stats(count*groups+loc2, :) = boxplot_stats_mmy(deltaT_tree(ts_s(count):ts_e(count),:,:),def_value)
      ;
      ; SM64cm_stats(count*groups+loc1, :) = boxplot_stats_mmy(SM64cm(ts_s(count):ts_e(count),:,:),def_value)
      ; SM64cm_stats(count*groups+loc2, :) = boxplot_stats_mmy(SM64cm_tree(ts_s(count):ts_e(count),:,:),def_value)
    end do ; count

    delete(Qle)
    ; delete(Qh)
    ; delete(deltaT)
    ; delete(Fwsoil)
    ; delete(SM64cm)

    delete(Qle_tree)
    ; delete(Qh_tree)
    ; delete(Fwsoil_tree)
    ; delete(deltaT_tree)
    ; delete(SM64cm_tree)
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

  Qle_stats   = where( ismissing(Qle_stats), -9999., Qle_stats)
  print(Qle_stats)
  Qle_stats@_FillValue = -9999.
  option@fout = "Qle_box_stats.txt"
  format      = case_sum+"f20.4"
  write_matrix (Qle_stats, format, option)
  delete(option@fout)

  ; Qh_stats   = where( ismissing(Qh_stats), -9999., Qh_stats)
  ; print(Qh_stats)
  ; Qh_stats@_FillValue = -9999.
  ; option@fout = "Qh_box_stats.txt"
  ; format      = case_sum+"f20.4"
  ; write_matrix (Qh_stats, format, option)
  ; delete(option@fout)

  ; Fwsoil_stats   = where( ismissing(Fwsoil_stats), -9999., Fwsoil_stats)
  ; print(Fwsoil_stats)
  ; Fwsoil_stats@_FillValue = -9999.
  ; option@fout = "Fwsoil_box_stats.txt"
  ; format      = case_sum+"f20.4"
  ; write_matrix (Fwsoil_stats, format, option)
  ; delete(option@fout)

  ; deltaT_stats   = where( ismissing(deltaT_stats), -9999., deltaT_stats)
  ; print(deltaT_stats)
  ; deltaT_stats@_FillValue = -9999.
  ; option@fout = "deltaT_box_stats.txt"
  ; format      = case_sum+"f20.4"
  ; write_matrix (deltaT_stats, format, option)
  ; delete(option@fout)

  ; SM64cm_stats   = where( ismissing(SM64cm_stats), -9999., SM64cm_stats)
  ; print(SM64cm_stats)
  ; SM64cm_stats@_FillValue = -9999.
  ; option@fout = "SM64cm_box_stats.txt"
  ; format      = case_sum+"f20.4"
  ; write_matrix (SM64cm_stats, format, option)
  ; delete(option@fout)

else if operation .eq. "plot" then

  txt_name  = "Qle_box_stats.txt"
  Qle_stats = asciiread(txt_name,(/total_box,5/),"float")

  txt_name  = "Qh_box_stats.txt"
  Qh_stats = asciiread(txt_name,(/total_box,5/),"float")

  txt_name  = "Fwsoil_box_stats.txt"
  Fwsoil_stats = asciiread(txt_name,(/total_box,5/),"float")

  txt_name  = "deltaT_box_stats.txt"
  deltaT_stats = asciiread(txt_name,(/total_box,5/),"float")

  txt_name  = "SM64cm_box_stats.txt"
  SM64cm_stats = asciiread(txt_name,(/total_box,5/),"float")

  ; ================================ Plotting ==================================
  x = ispan(1,total_box,1)

  wks = gsn_open_wks("pdf",pic_name)   ;"./plots/"+            ; send graphics to PNG file

  res            = True                         ; plot mods desired
  res@tmXBLabels = (/"All","","","","", "2001-09","","","","",\
                     "2017-19","","","","","2010-11","","",""/) ; labels for each box
  ; resources for polylines that draws the boxes
  llres                   = True
  llres@gsLineThicknessF  = 3.                ; line thickness

  ; resources that control color and width of boxes
  opti           = True
  opti@boxWidth  = 0.9
  opti@boxColors = (/"blue","red","green","orange","black", "blue","red","green","orange","black",\
                     "blue","red","green","orange","black", "blue","red","green","orange"/)

  ; plot = boxplot_mmy(wks,x,Qle_stats,opti,res,llres, 0, 150, "Qle")
  ; draw(wks)                                    ; boxplot does not call these
  ; frame(wks)                                   ; for you
  ; plot = boxplot_mmy(wks,x,Qh_stats,opti,res,llres, 0, 250, "Qh")
  ; draw(wks)                                    ; boxplot does not call these
  ; frame(wks)                                   ; for you
  plot = boxplot_mmy(wks,x,Fwsoil_stats,opti,res,llres, 0, 1.1, "Fwsoil")
  draw(wks)                                    ; boxplot does not call these
  frame(wks)                                   ; for you
  ; plot = boxplot_mmy(wks,x,deltaT_stats,opti,res,llres, 0, 6, "deltaT")
  ; draw(wks)                                    ; boxplot does not call these
  ; frame(wks)                                   ; for you
  ; ; plot = boxplot_mmy(wks,x,WatTable_stats,opti,res,llres, 0, 15, "WatTable")
  ; ; draw(wks)                                    ; boxplot does not call these
  ; ; frame(wks)                                   ; for you
  ; plot = boxplot_mmy(wks,x,SM64cm_stats,opti,res,llres, 0.1, 0.4, "SM64cm")
  ; draw(wks)                                    ; boxplot does not call these
  ; frame(wks)                                   ; for you

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
end if ; operation
end if ; operation


end
