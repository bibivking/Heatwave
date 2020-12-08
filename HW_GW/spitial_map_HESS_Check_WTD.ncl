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
year        = 2009
filename_hw = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL-30x3+30yr/outputs/cable_out_"+year+"_SE_Aus.nc"
f           = addfile (filename_hw,"r")
var         = dim_avg_n_Wrap(f->WatTable(:,:,:),0)
delete(f)
filename_hw = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL-30x2+30yr/outputs/cable_out_"+year+"_SE_Aus.nc"
f           = addfile (filename_hw,"r")
var2        = dim_avg_n_Wrap(f->WatTable(:,:,:),0)
delete(f)
; ================================= Plotting ===================================
res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                    ;��ɫ
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

pic1 = "./plots/Check_WTD_3_"+year
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"GMT_haxby");"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

levels = fspan(0,20,21)
res@cnLevels = levels
res@txString = "Water Table Depth (CTL-30+30yr)" ;var@long_name
plots = gsn_csm_contour(wks1,var,res)
res@txString = "Water Table Depth (CTL-30x2+30yr)" ;var@long_name
plots = gsn_csm_contour(wks1,var2,res)
; delete(levels)
; delete(res@cnLevels)
; levels = fspan(-10,10,21)
; res@cnLevels = levels
; plots = gsn_csm_contour(wks1,var-var2,res)
end
