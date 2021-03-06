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

; ================================== Options ===================================

filename_hw = "/g/data/w35/Shared_data/Observations/Fluxnet_data/Post-processed_PLUMBER2_outputs/Nc_files/Flux/AU-Tum_2002-2017_OzFlux_Flux.nc"
f           = addfile (filename_hw,"r")
var         = f->Qle(:,0,0) ;*24*2600
var1        = f->Qh(:,0,0) ;*24*2600
var2        = f->Rnet(:,0,0) ;*24*2600
delete(f)

; var2 = new((/12,12,14/),float)
; do month = 1,12
;   if month .lt. 10 then
;     filename_hw = "/g/data/w35/mm3972/data/CERES/CERES_EBAF_Ed4.1_Subset_20130"+month+"-"+"20130"+month+".nc"
;   else
;     filename_hw = "/g/data/w35/mm3972/data/CERES/CERES_EBAF_Ed4.1_Subset_2013"+month+"-"+"2013"+month+".nc"
;   end if
;
;   f           = addfile (filename_hw,"r")
;   print(f->sfc_net_tot_clr_t_mon(0,{-40:-28},{140:154}))
;   var2(month-1,:,:) = f->sfc_net_tot_clr_t_mon(0,{-40:-28},{140:154})
;   delete(f)
; end do

; filename_hw = "/g/data/w35/Shared_data/AWAP_3h_v1/SWdown/AWAP.SWdown.3hr.2013.nc"
; f           = addfile (filename_hw,"r")
; var2 = f->SWdown(:,{-40:-28},{140:154})
; delete(f)


; filename_hw = "/g/data/w35/Shared_data/Observations/AWAP_all_variables/daily/rad/AWAP_daily_rad_1990_2019.nc"
; f           = addfile (filename_hw,"r")
; var2        = f->rad(:,{-40:-28},{140:154})
; delete(f)
; ================================= Plotting ===================================
res                    = True
; res@cnFillMode         = "RasterFill"            ; Raster Mode
; res@cnFillOn           = True                    ;��ɫ
; res@tmBorderThicknessF = 3.0

res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ����ҳ
res@lbLabelBarOn       = True
;
; ;************** ����labelbar ***************
; res@lbBoxLinesOn       = True                       ;�ر�lbar box ����
; res@lbTitleFont        = 25
; res@lbLabelFont        = 25
; res@lbTitleFontHeightF = 0.013
; res@lbLabelFontHeightF = 0.013
; res@txString           = ""
; res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
; res@tmYLLabelFont      = 25
;
; ;*************** ���õ�ֵ�� *****************
; res@cnLinesOn          = False                       ; �رյ�ֵ������
; res@cnLineColor        = "black"
; res@cnLineThicknessF   = 1.5
; res@cnLineLabelsOn     = False
; res@gsnMaximize        = True
; res@cnExplicitLabelBarLabelsOn = True   ;?
; res@cnLevelSelectionMode = "ExplicitLevels"
;
; ;************ ����������ֵ��Χ **************
; res@tmXBLabelFontThicknessF = 0.015
; res@tmYLLabelFontThicknessF = 0.015
; res@tmXBLabelFontHeightF = 0.015
; res@tmYLLabelFontHeightF = 0.015
; res@tmYLMode  = "Explicit"
; res@tmXBMode  = "Explicit"
;
; res@gsnLeftString       = ""
; res@gsnRightString      = ""

pic1 = "./plots/Check_Tum-Qle-Qh"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"srip_reanalysis");"WhiteBlueGreenYellowRed");"GMT_haxby");"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

x = ispan(0,23,1) ;"Water Table Depth" ;var@long_name
; levels = fspan(-1.,1.,21)
; res@cnLevels = levels
tmp  = var(0:23)
tmp1 = var1(0:23)
tmp2 = var2(0:23)
do counter = 1,30
  tmp  = tmp + var(counter*24:counter*24+23)
  tmp1 = tmp1 + var1(counter*24:counter*24+23)
  tmp2 = tmp2 + var2(counter*24:counter*24+23)
end do

tmp  = tmp/31
tmp1 = tmp1/31
tmp2 = tmp2/31

res@txString = "Qle"
plots = gsn_csm_xy(wks1,x,tmp,res)
res@txString = "Rnet"
plots = gsn_csm_xy(wks1,x,tmp1,res)
res@txString = "Qh"
plots = gsn_csm_xy(wks1,x,tmp2,res)
  ; plots = gsn_csm_contour(wks1,var2(day_num, :,:),res)
; end do
; counter = 0
; do year = 1990, 2019
;   if mod(year,4) .eq. 0 then
;     day_sum = 366
;   else
;     day_sum = 365
;   end if
;   do day_num = 0, day_sum -1
;     res@gsnLeftString = year+"-"+day_num
;     if day_num .lt. 60 .or. day_num .gt. day_sum - 61 then
;       plots = gsn_csm_contour(wks1,var2(counter,:,:),res)
;       counter = counter + 1
;     else
;       counter = counter + 1
;     end if
;   end do
; end do

; delete(levels)
; delete(res@cnLevels)
; levels = fspan(-10,10,21)
; res@cnLevels = levels
; plots = gsn_csm_contour(wks1,var-var2,res)
end
