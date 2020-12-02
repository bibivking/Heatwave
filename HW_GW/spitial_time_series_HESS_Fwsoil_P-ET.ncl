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

; ================================= Predefine ==================================
pic_name  = "Fwsoil_P-ET"
case_name = (/"CTL","FREE_DRAIN"/) ;,"PUMPx10","PUMPx100"/);,"PUMP_IRRG","PUMP_IRRGx10","PUMP_IRRGx100"/)

PFT       = "all" ; "tree" ;"all"
plot_type = "all" ; "spitial" ; "time_series"
scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
case_sum  = dimsizes(case_name)

year_s   = 2000
year_e   = 2019
year_sum = year_e-year_s+1
total    = year_sum*12
plot_sum = case_sum*2

; set PFT
if PFT .eq. "all" then
  iveg_s = 0
  iveg_e = 20
else if PFT .eq. "tree" then
  iveg_s = 1
  iveg_e = 4
end if
end if

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

; =================================== Path =====================================
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"+case_name(case_num)+"/outputs-spinup200yr"
end do

; ================================== SET Dim ===================================
filename  = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f         = addfile (filename,"r")

lat_size  = dimsizes(f->y)
lon_size  = dimsizes(f->x)
CABLE_lat = f->y
CABLE_lon = f->x
default_value = f->Rainf@_FillValue
lat_sum = dimsizes(f->Rainf(0,{A:B},0))
lon_sum = dimsizes(f->Rainf(0,0,{C:D}))

Fwsoil   = new((/case_sum, total, lat_sum, lon_sum/), float, default_value)
Rain     = new((/case_sum, total, lat_sum, lon_sum/), float, default_value)
Evap     = new((/case_sum, total, lat_sum, lon_sum/), float, default_value)
Fwsoil   = 0.
Rain     = 0.
Evap     = 0.

; PFT mask
PFT_mask  = where((f->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
PFT_mask_multi = new((/total,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, total-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do

delete(PFT_mask)
delete(f)

; ================================= Read Data ==================================
do case_num = 0,case_sum-1
  ts = 0
  do year = year_s,year_e
    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if
    filename1  = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
    f1         = addfile (filename1,"r")
    day_s      = 0
    Fwsoil_tmp = f1->Fwsoil(:,{A:B},{C:D})
    Rain_tmp   = f1->Rainf(:,{A:B},{C:D})
    Evap_tmp   = f1->Evap(:,{A:B},{C:D})
    do month = 0, 11
      Fwsoil(case_num,ts,:,:) = dim_avg_n_Wrap(Fwsoil_tmp(day_s:day_s+dom(month)-1,:,:),0)
      Rain(case_num,ts,:,:)   = dim_sum_n_Wrap(Rain_tmp(day_s:day_s+dom(month)-1,:,:),0)*3600*24
      Evap(case_num,ts,:,:)   = dim_sum_n_Wrap(Evap_tmp(day_s:day_s+dom(month)-1,:,:),0)*3600*24
      day_s                   = day_s + dom(month)
      ts                      = ts + 1
    end do ; month
    delete(Fwsoil_tmp)
    delete(Rain_tmp)
    delete(Evap_tmp)
    delete(f1)
  end do ; year
end do ; case_num
delete(ts)

do case_num = 0, case_sum - 1
    Fwsoil(case_num,:,:,:) = where(PFT_mask_multi .eq. 1, Fwsoil(case_num,:,:,:), default_value)
    Rain(case_num,:,:,:)   = where(PFT_mask_multi .eq. 1, Rain(case_num,:,:,:), default_value)
    Evap(case_num,:,:,:)   = where(PFT_mask_multi .eq. 1, Evap(case_num,:,:,:), default_value)
end do
delete(PFT_mask_multi)

; ============================= plot spitial map ===============================
if (plot_type .eq. "spitial" ) .or. (plot_type .eq. "all") then

    Fwsoil_spitial = dim_avg_n_Wrap(Fwsoil,1)
    Rain_spitial   = dim_sum_n_Wrap(Rain,1)/year_sum
    Evap_spitial   = dim_sum_n_Wrap(Evap,1)/year_sum

    ; ______________________________ Plotting ______________________________
    pic1 = "./plots/HESS_spitial_"+plot_type+"_PFT-"+PFT+"_"+year_s+"-"+year_e
    wks1 = gsn_open_wks("pdf",pic1)
    gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    res1                    = True
    res1@cnFillMode         = "RasterFill"            ; Raster Mode
    res1@cnFillOn           = True                            ;��ɫ
    res1@tmBorderThicknessF = 3.0

    res1@gsnDraw            = False  ; Don't draw plots
    res1@gsnFrame           = False  ; ����ҳ
    res1@lbLabelBarOn       = True

    ;************** ����labelbar ***************
    res1@lbBoxLinesOn       = True                       ;�ر�lbar box ����
    res1@lbTitleFont        = 25
    res1@lbLabelFont        = 25
    res1@lbTitleFontHeightF = 0.013
    res1@lbLabelFontHeightF = 0.013
    res1@txString           = ""
    res1@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
    res1@tmYLLabelFont      = 25

    ;*************** ���õ�ֵ�� *****************
    res1@cnLinesOn          = False                       ; �رյ�ֵ������
    res1@cnLineColor        = "black"
    res1@cnLineThicknessF   = 1.5
    res1@cnLineLabelsOn     = False
    res1@gsnMaximize        = True
    res1@cnExplicitLabelBarLabelsOn = True   ;?
    res1@cnLevelSelectionMode = "ExplicitLevels"

    ;************ ����������ֵ��Χ **************
    res1@tmXBLabelFontThicknessF = 0.015
    res1@tmYLLabelFontThicknessF = 0.015
    res1@tmXBLabelFontHeightF = 0.015
    res1@tmYLLabelFontHeightF = 0.015
    res1@tmYLMode  = "Explicit"
    res1@tmXBMode  = "Explicit"
    res1@gsnLeftString = ""
    res1@gsnRightString = ""

    ;**************** ͼ�������� *****************
    pres1                    = True                                      ; ���� panel res1ources.
    pres1@gsnMaximize        = True
    pres1@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
    pres1@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
    pres1@lbBoxLinesOn       = True                       ;�ر�lbar box ����
    pres1@lbTitleFont        = 25
    pres1@lbLabelFont        = 25
    ;pres1@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
    pres1@lbTitleFontHeightF = 0.01
    pres1@lbLabelFontHeightF = 0.01

    plots = new(plot_sum,graphic)

    ; Fwsoil
    res1@cnFillPalette = "WhiteBlueGreenYellowRed"
    res1@cnLevels      = fspan(0,1,21)
    do case_num = 0,case_sum-1
      res1@tiMainString = case_name(case_num) + " Fwsoil (-)"
      plots(case_num) = gsn_csm_contour(wks1,Fwsoil_spitial(case_num,:,:)*1.,res1)
    end do

    ; P - ET
    res1@cnFillPalette = "ViBlGrWhYeOrRe"
    res1@cnLevels = fspan(-200,200,21)
    do case_num = 0,case_sum-1
      res1@tiMainString = case_name(case_num) + " P - ET (mm/yr)"
      plots(case_sum + case_num) = gsn_csm_contour(wks1,Rain_spitial(case_num,:,:)\
                                                   -Evap_spitial(case_num,:,:),res1)
    end do
    gsn_panel(wks1,(/plots/),(/1,plot_sum/),pres1)
    delete(Fwsoil_spitial)
    delete(Rain_spitial)
    delete(Evap_spitial)
    delete(plots)
end if

; ============================= plot time series ===============================
if (plot_type .eq. "time_series") .or. (plot_type .eq. "all") then

    Fwsoil_time_series = dim_avg_n_Wrap(Fwsoil,(/2,3/))
    Rain_time_series   = dim_avg_n_Wrap(Rain,(/2,3/))
    Evap_time_series   = dim_avg_n_Wrap(Evap,(/2,3/))

    Rain_acl_time_series = dim_cumsum_n_Wrap(Rain_time_series, 2, 1)
    Evap_acl_time_series = dim_cumsum_n_Wrap(Evap_time_series, 2, 1)
    ; opt = 2 : treat _FillValue as numeric zero when summing.

    Fwsoil_t_s_sep   = new((/case_sum*6,total/),float, default_value)
    P_ET_t_s_sep     = new((/case_sum*6,total/),float, default_value)
    P_ET_acl_t_s_sep = new((/case_sum*6,total/),float, default_value)

    Fwsoil_t_s_sep   = default_value
    P_ET_t_s_sep     = default_value
    P_ET_acl_t_s_sep = default_value

    case_counter = 0
    do spell_num = 0, 5
      if spell_num .eq. 0 then ; 2000 - 2000
        month_s = (2000 - year_s)*12
        month_e = (2000 - year_s)*12 + 12
      else if spell_num .eq. 1 then ; 2001 - 2009
        month_s = (2001 - year_s)*12
        month_e = (2009 - year_s)*12 + 12
      else if spell_num .eq. 2 then ; 2010 - 2011
        month_s = (2010 - year_s)*12
        month_e = (2011 - year_s)*12 + 12
      else if spell_num .eq. 3 then ; 2012 - 2015
        month_s = (2012 - year_s)*12
        month_e = (2015 - year_s)*12 + 12
      else if spell_num .eq. 4 then ; 2016 - 2016
        month_s = (2016 - year_s)*12
        month_e = (2016 - year_s)*12 + 12
      else if spell_num .eq. 5 then ; 2017 - 2019
        month_s = (2017 - year_s)*12
        month_e = (2019 - year_s)*12 + 11
      end if
      end if
      end if
      end if
      end if
      end if
      do case_num = 0, case_sum-1
        Fwsoil_t_s_sep(case_counter, month_s:month_e) = Fwsoil_time_series(case_num, month_s:month_e)
        P_ET_t_s_sep(case_counter, month_s:month_e)   = Rain_time_series(case_num, month_s:month_e) \
                                                      - Evap_time_series(case_num, month_s:month_e)
        P_ET_acl_t_s_sep(case_counter, month_s:month_e) = Rain_acl_time_series(case_num, month_s:month_e) \
                                                        - Evap_acl_time_series(case_num, month_s:month_e)
        case_counter = case_counter + 1
      end do ; case_num
    end do ; spell_num

    ; ______________________________ Plotting ______________________________
    pic2 = "./plots/HESS_time_series_"+plot_type+"_PFT-"+PFT+"_"+year_s+"-"+year_e
    wks2 = gsn_open_wks("pdf",pic2)
    gsn_define_colormap(wks2,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    res2                      = True
    res2@gsnMaximize          = True
    res2@gsnDraw              = False  ; Don't draw plots
    res2@gsnFrame             = False  ; ����ҳ
    res2@vpHeightF            = 0.3
    res2@vpWidthF             = 0.90
    res2@vpXF                 = 0.08
    res2@vpYF                 = 0.98
    res2@pmLegendWidthF       = 0.2                     ; res2ize legend width
    res2@pmLegendHeightF      = 0.2                     ; and height
    res2@tiXAxisFontHeightF   = 0.02
    res2@tiYAxisFontHeightF   = 0.02
    res2@tmXBLabelFontHeightF = 0.02
    res2@tmYLLabelFontHeightF = 0.02
    res2@gsnPaperOrientation  = "portrait";"portrait";"landscape" ;"auto"

    res2@tmXBMode             = "Explicit"
    res2@tmBorderThicknessF   = 3.0
    res2@tmXBValues           = ispan(1,total,12)
    res2@tmXBLabels           = tostring(ispan(year_s,year_e,1))
    res2@pmLegendDisplayMode  = "never" ;"Always"

    res2@xyLineThicknesses    = (/ 1., 1., 1., 1., 1., 1., \
                                   1., 1., 1., 1., 1., 1. /)
    ; res2@xyDashPatterns       = (/ 0, 5, 0, 5, 0, 5, 0, 5, 0, 5, 0, 5 /)
    res2@xyDashPatterns       = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
    res2@xyLineColors         = (/ "black","gray", "red","orange", "blue", "green", \
                                   "black","gray", "blue", "green", "red","orange" /)
    ; res2@xyExplicitLegendLabels = case_name

    res2@tiXAxisString        = "months"
    ;**************** ͼ�������� *****************
    pres2                    = True                                      ; ���� panel res2ources.
    pres2@gsnMaximize        = True
    pres2@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
    pres2@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
    pres2@lbBoxLinesOn       = True                       ;�ر�lbar box ����
    pres2@lbTitleFont        = 25
    pres2@lbLabelFont        = 25
    pres2@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
    pres2@lbTitleFontHeightF = 0.2
    pres2@lbLabelFontHeightF = 0.2

    x             = ispan(1,total,1)
    plots         = new(3,graphic)

    res2@tiYAxisString = "Fwsoil (-)"
    plots(0)  = gsn_csm_xy(wks2,x,Fwsoil_t_s_sep,res2)
    res2@tiYAxisString = "P - ET (mm/mth)"
    plots(1)  = gsn_csm_xy(wks2,x,P_ET_t_s_sep,res2)
    res2@tiYAxisString = "Accumlation of P - ET (mm)"
    plots(2)  = gsn_csm_xy(wks2,x,P_ET_acl_t_s_sep,res2)

    gsn_panel(wks2,(/plots/),(/3,1/),pres2)
    delete(plots)
end if

end
