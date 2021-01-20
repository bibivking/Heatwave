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

undef("nd_to_1d_mmy")
function nd_to_1d_mmy(var[*][*][*]:numeric,def_value:numeric)
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

  return(var_masked)
end

; ================================== Options ===================================
year_s       = 2000
year_e       = 2019

scale        = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
heatwave     = "hw-90th" ; "hw-90th" ; "all-day"
case_name    = (/"CTL-30x3+30yr/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr"/) ;, "FREE_DRAIN/outputs-spinup30x3yr", \
                 ;"CTL-30x3+30yr-DeepRoot/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
case_name_brf= "CTL" ;_vs_FD_vs_DR"
var_name     = "Dr";"Qrec";"deltaT"
operation    = "calc"    ;"calc" ; "plot"

order_name   = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

iveg_tree_s    = 1
iveg_tree_e    = 4
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

; _____________________________ Calc Bin Stats ________________________________
;  all years,  2001-2009,  2010-2011,   2017-2019
;   2000-1-1,   2001-1-1,   2010-1-1,    2017-1-1
ts_s    = (/        0,         366,       3653,       6210/)
; 2019-12-31, 2009-12-31, 2011-12-31,  2019-12-31
ts_e    = (/day_sum-1,        3652,       4382,       7304/)
periods = (/ "2000-19",  "2001-09",  "2010-11",  "2017-19"/)


if operation .eq. "calc" then
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

  do case_num = 0,0;case_sum-1
    print("case_num = " +case_num)
    ; ______________________________ SET DATA _________________________________
    WTD      = new((/day_sum,lat_sum,lon_sum/),float,def_value)
    WTD      = 0.
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
      filename2 = path(case_num+1)+"/cable_out_"+year+"_SE_Aus.nc"
      f2 = addfile (filename2,"r")


      WTD(day_step_s:day_step_e,:,:)   = f1->WatTable(:,{A:B},{C:D})

      if var_name .eq. "Qle" then
        var(day_step_s:day_step_e,:,:)   = f1->$var_name$(:,{A:B},{C:D})
      else if var_name .eq. "Qh" then
        var(day_step_s:day_step_e,:,:)   = f1->$var_name$(:,{A:B},{C:D})
      else if var_name .eq. "Fwsoil" then
        var(day_step_s:day_step_e,:,:)   = f1->$var_name$(:,{A:B},{C:D})
      else if var_name .eq. "deltaT" then
        var(day_step_s:day_step_e,:,:)   = f1->VegT(:,{A:B},{C:D}) - f1->Tair(:,{A:B},{C:D})
      else if var_name .eq. "SM64cm" then
        var(day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022  \
                                            + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
                                            + f1->SoilMoist(:,2,{A:B},{C:D})*0.154 \
                                            + f1->SoilMoist(:,3,{A:B},{C:D})*0.409)\
                                            / 0.643
      else if var_name .eq. "Dr" then
        Qrec_tmp  = f1->Qrecharge(:,{A:B},{C:D})
        Qrec_tmp2 = f2->Qrecharge(:,{A:B},{C:D})
        Qrec_tmp =  where(Qrec_tmp  .ge. 0 .and. Qrec_tmp  .ne. def_value, Qrec_tmp*24*3600, def_value)
        Qrec_tmp2 = where(Qrec_tmp2 .ge. 0 .and. Qrec_tmp2 .ne. def_value, Qrec_tmp2*24*3600, def_value)
        printVarSummary(Qrec_tmp)
        printVarSummary(Qrec_tmp2)
        var(day_step_s:day_step_e,:,:) = Qrec_tmp2 - Qrec_tmp
        delete(Qrec_tmp)
        delete(Qrec_tmp2)
        WTD(day_step_s:day_step_e,:,:) = where(.not. ismissing(var(day_step_s:day_step_e,:,:)),\
                                              WTD(day_step_s:day_step_e,:,:), def_value)
      else if var_name .eq. "Qrec" then
        var(day_step_s:day_step_e,:,:) = f1->Qrecharge(:,{A:B},{C:D})*24*3600
      end if
      end if
      end if
      end if
      end if
      end if
      end if

      ; calc fist day of next year
      day_step_s = day_step_e + 1

      delete(filename1)
      delete(f1)
      delete(filename2)
      delete(f2)
    end do ; year

    ; _____________________________ MASK HW DAY ________________________________
    if heatwave .eq. "hw-90th" then
      print("mask hw-90th")
      WTD      = where(HW .eq. 1, WTD, def_value)
      var      = where(HW .eq. 1, var, def_value)
    end if

    ; _______________________________ MASK PFT __________________________________
    print("mask PFT")
    ; tree
    WTD_tree = where(PFT_mask_multi .eq. 1, WTD, def_value)
    var_tree = where(PFT_mask_multi .eq. 1, var, def_value)
    ; non-trees
    WTD      = where(PFT_mask_multi .eq. 2, WTD, def_value)
    var      = where(PFT_mask_multi .eq. 2, var, def_value)

    loc1 = case_num
    loc2 = case_num + 1 ;case_sum

    do count = 0, dimsizes(ts_s)-1 ; four periods

      ; ______________________________ Trees _______________________________
      WTD_tmp = nd_to_1d_mmy(WTD_tree(ts_s(count):ts_e(count),:,:),def_value)
      tmp     = nd_to_1d_mmy(var_tree(ts_s(count):ts_e(count),:,:),def_value)

      print("write out")
      option        = True
      option@row    = False

      var_tmp     = new((/dimsizes(WTD_tmp), 2/),float,def_value)
      var_tmp(:,0)= WTD_tmp
      var_tmp(:,1)= tmp
      option@fout = "scatter_WTD_"+var_name+"_"+case_name_brf+"_PFT-tree_"+periods(count)+".txt"
      format      = "2f20.4"
      write_matrix (var_tmp, format, option)
      delete(option)
      delete(var_tmp)
      delete(WTD_tmp)
      delete(tmp)

      ; ______________________________ Nontrees _______________________________
      WTD_tmp = nd_to_1d_mmy(WTD(ts_s(count):ts_e(count),:,:),def_value)
      tmp     = nd_to_1d_mmy(var(ts_s(count):ts_e(count),:,:),def_value)

      option        = True
      option@row    = False
      var_tmp     = new((/dimsizes(WTD_tmp), 2/),float,def_value)
      var_tmp(:,0)= WTD_tmp
      var_tmp(:,1)= tmp
      option@fout = "scatter_WTD_"+var_name+"_"+case_name_brf+"_PFT-nontree_"+periods(count)+".txt"
      format      = "2f20.4"
      write_matrix (var_tmp, format, option)
      delete(option)
      delete(var_tmp)
      delete(WTD_tmp)
      delete(tmp)

    end do ; count

    delete(WTD)
    delete(WTD_tree)
    delete(var)
    delete(var_tree)
  end do ; case_num

  if heatwave .eq. "hw-90th" then
    delete(HW)
  end if
  delete(PFT_mask_multi)

else if operation .eq. "plot" then

  ; ============================= setting plot ===============================
  res                        = True
  res@gsnDraw                = True  ; Don't draw plots
  res@gsnFrame               = True  ; ����ҳ
  res@gsnMaximize            = True
  res@gsnPaperOrientation = "portrait"

  res@vpHeightF              = 0.90
  res@vpWidthF               = 0.90
  res@vpXF                   = 0.08
  res@vpYF                   = 0.98

  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Right";"Top"               ; Change location of
  res@pmLegendParallelPosF   = .25                 ; move units right
  res@pmLegendOrthogonalPosF = -0.4                ; move units down
  res@pmLegendWidthF         = 0.125               ; Change width and
  res@pmLegendHeightF        = 0.15                ; height of legend.
  res@lgPerimOn              = False                ; turn off/on box around
  res@lgLabelFontHeightF     = .015                ; label font height
  res@tiYAxisString          = "~F5~D"+"~F25~T (~S~o~N~C)"

  ; res@tmXBBorderOn      = False
  ; res@tmXBOn            = False
  ; res@tmXUseBottom      = False
  ; res@tmXTLabelsOn      = True
  ; res@tmYLMinorPerMajor = 1
  ; res@tmYRBorderOn      = False
  ; res@tmYROn            = False
  ; res@tmYUseLeft        = False
  ; res@trYReverse        = True

  res@trXMinF             = 0.
  res@trXMaxF             = 15.

  res@xyMarkLineMode      = "Markers"
  res@xyMarker            = 16                     ; choose type of marker
  res@xyMarkerOpacityF    = 0.01
  res@xyMarkerColor       = "lightskyblue";"navy"
  res@xyMarkerSizeF       = 0.01
  res@xyExplicitLegendLabels = "GW"  ; create explicit labels

  if var_name .eq. "Qle" then
    res@trYMinF = 0.
    res@trYMaxF = 300.
  else if var_name .eq. "Qh" then
    res@trYMinF = -50.
    res@trYMaxF = 350.
  else if var_name .eq. "Fwsoil" then
    res@trYMinF = 0.
    res@trYMaxF = 1.
  else if var_name .eq. "deltaT" then
    print("skip")
    ; res@trYMinF = -3.
    ; res@trYMaxF = 6.
  else if var_name .eq. "SM64cm" then
    res@trYMinF = 0.05
    res@trYMaxF = 0.45
  end if
  end if
  end if
  end if
  end if

  txres                   = True
  txres@txPerimOn         = False
  txres@txFont            = 25
  txres@txFontHeightF     = 0.03

  amres                   = True
  amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
  amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
  amres@amJust            = "TopLeft"


  ; ================================ Plot =====================================
  do count = 1, dimsizes(ts_s)-1

    pic_name     = "HW_scatter_"+var_name+"_WTD_"+"_"+heatwave+"_"+scale+"_"+case_name_brf+periods(count)
    wks = gsn_open_wks("pdf",pic_name)
    gsn_define_colormap(wks,"rainbow")

    ; =============================== READ IN ======================================
    txt_name = "./txt/scatter_WTD_"+var_name+"_"+case_name_brf+"_PFT-tree_"+periods(count)+".txt"
    data_tmp = asciiread(txt_name, -1,"float") ; read all data as one-D
    print(data_tmp(0:3))
    data     = onedtond(data_tmp, (/dimsizes(data_tmp)/2,2/))
    print(data(0,:))
    res@tiXAxisString = "WTD (m)"
    plot    = gsn_csm_xy (wks, data(:,0), data(:,1), res)
    txid1   = gsn_create_text(wks, order_name(count), txres)
    annoid1 = gsn_add_annotation(plot, txid1, amres)
    delete(data_tmp)
    delete(data)
    delete(pic_name)
    delete(wks)
  end do

end if
end if

exit
end
