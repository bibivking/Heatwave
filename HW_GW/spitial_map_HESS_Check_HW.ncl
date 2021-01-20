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
year_s    = 2019
year_e    = 2019

filename_hw = "/g/data/w35/mm3972/scripts/ehfheatwaves/HW_Event_Indicator_1970-2019.nc"
;"/g/data/w35/mm3972/scripts/ehfheatwaves/EHF_heatwaves____daily.nc"
;"/g/data/w35/mm3972/scripts/ehfheatwaves/EHI_heatwaves____daily.nc"
f    = addfile (filename_hw,"r")

; tmp  = doubletofloat(f->event(0,:,:))*0
; do year = year_s, year_e
;   time_s    = (year-1970)*365
;   time_e    = (year-1970+1)*365-1
;   tmp       = tmp + dim_sum_n_Wrap(doubletofloat(f->event(time_s:time_s+58,:,:)),0) \
;                   + dim_sum_n_Wrap(doubletofloat(f->event(time_e - 30 :time_e,:,:)),0)
; end do

; tmp   = f->event(0,:,:)
do year = year_s, year_e

  time_s    = (year-1970)*365 + 12
  time_e    = (year-1970+1)*365-1 + 12
  tmp   = f->event(time_s:time_e,:,:)
  ; tmp       = tmp + dim_sum_n_Wrap(f->event(time_s:time_s+58,:,:),0) \
  ;                 + dim_sum_n_Wrap(f->event(time_e - 30 :time_e,:,:),0)

end do

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
res@txString      = ""
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
res@tmYLMode  = "Explicit"
res@tmXBMode  = "Explicit"

res@gsnMaximize         = True
res@gsnLeftString       = ""
res@gsnRightString      = ""

pic1 = "./plots/Check_heatwave"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

res@txString = "total Heatwave days"
levels = fspan(0,1,21)
res@cnLevels = levels
do time = 0, 60
  plots = gsn_csm_contour(wks1,tmp(time,:,:),res)
end do


end
