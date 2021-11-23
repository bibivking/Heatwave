
;MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
; Author  : Mengyuan Mu
; Version : 1.0 (22.03.2021)"
; Email   : mu.mengyuan815@gmail.com
;WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


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
; ; 2009
; hw_event  = "hw2009_15Oct"
; year_s    = 2008
; year_e    = 2009
; doy_s     = 335 ; 1 Dec
; doy_e     = 58  ; 28 Feb

; ; 2013
; hw_event  = "hw2013_15Oct"
; year_s    = 2012
; year_e    = 2013
; doy_s     = 335 ; 1 Dec
; doy_e     = 30  ; 31 Jan

; 2019
hw_event  = "hw2019_15Oct"
year_s    = 2018
year_e    = 2019
doy_s     = 334 ; 1 Dec
doy_e     = 58  ; 28 Feb

; ===================================== Path ====================================
filename_hw = "/g/data/w35/mm3972/scripts/ehfheatwaves/nc_file/AWAP_daily_tmean_1970_2019.nc"
;"/g/data/w35/Shared_data/Observations/AWAP_all_variables/daily/tmax/AWAP_daily_tmax_1970_2019.nc"
;"/g/data/w35/mm3972/scripts/ehfheatwaves/nc_file/SE_AUS/HW_Event_Indicator_1970-2019.nc"
f    = addfile (filename_hw,"r")

; time_s = doy_s
; time_e = doy_e

; do year = 1970, year_s-1
;   if ( mod(year,4) .eq. 0) then
;       time_s = time_s + 366
;   else
;       time_s = time_s + 365
;   end if
; end do

; do year = 1970, year_e-1
;   if ( mod(year,4) .eq. 0) then
;       time_e = time_e + 366
;   else
;       time_e = time_e + 365
;   end if
; end do

; tmp = f->tmax(time_s:time_e,:,:)
; ;f->event(time_s:time_e,:,:)
; ;f->tmax(time_s:time_e,:,:)
; print(tmp&time)

; ================================= Plotting ===================================
res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0

res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ����ҳ
res@lbLabelBarOn       = True

;************** ����labelbar ***************
res@lbBoxLinesOn       = True                       ;�ر�lbar box ����
res@lbTitleFont        = 25
res@lbLabelFont        = 25
res@lbTitleFontHeightF = 0.013
res@lbLabelFontHeightF = 0.013
res@txString           = ""
res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont      = 25

;*************** ���õ�ֵ�� *****************
res@cnLinesOn          = False                       ; �رյ�ֵ������
res@cnLineColor        = "black"
res@cnLineThicknessF   = 1.5
res@cnLineLabelsOn     = False
res@gsnMaximize        = True
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode = "ExplicitLevels"

;************ ����������ֵ��Χ **************
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015
res@tmYLMode             = "Explicit"
res@tmXBMode             = "Explicit"

res@gsnMaximize         = True
res@gsnLeftString       = ""
res@gsnRightString      = ""

; pic = "./plots/19Oct/Heatwave_temp_" + hw_event
pic = "check" 
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")
plots = gsn_csm_contour(wks,f->tmean(0,:,:)*1.,res)

; levels = fspan(20,45,26)
; res@cnLevels = levels

; if mod(year_s,4) .eq. 0 then
;   leap_doy = 366
; else
;   leap_doy = 365
; end if

; counter = 0

; do time = 0, dimsizes(tmp(:,0,0))-1
;   print(time)
;   if ((doy_s + time) .gt. leap_doy-1) then
;     counter = counter + 1
;     res@gsnLeftString = year_e+" "+ counter
;     print(res@gsnLeftString)
;   else
;     res@gsnLeftString = year_s+" "+ (doy_s + time+1)
;     print(res@gsnLeftString)
;   end if
;   plots = gsn_csm_contour(wks,tmp(time,:,:),res)
; end do

end
