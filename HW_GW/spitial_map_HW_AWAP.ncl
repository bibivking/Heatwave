;*************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*************************************************************

begin

month      = 2
dom        = 28

pic = "AWAP_Tmax_Feb"

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"NCV_bright")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")
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

res@gsnAddCyclic       = False
res@trYMinF            = -44
res@trYMaxF            = -10
res@mpMinLatF          = -44
res@mpMaxLatF          = -10
res@trXMinF            = 112
res@trXMaxF            = 154
res@mpMinLonF          = 112
res@mpMaxLonF          = 154

res@tmYLMode  = "Explicit"             ;��Y����������������,��ʾ�Լ���Ҫ��ֵ
res@tmYLValues=(/-10,-20,-30,-40/)
res@tmYLLabels=(/"10~S~o~N~S","20~S~o~N~S","30~S~o~N~S","40~S~o~N~S"/)
res@tmXBMode  = "Explicit"
res@tmXBValues=(/120,130,140,150/)
res@tmXBLabels=(/"120~S~o~N~E","130~S~o~N~E","140~S~o~N~E","150~S~o~N~E"/)
res@cnLevels  = ispan(0,50,5)


ncname = "/g/data/w35/mm3972/data/AWAP/AWAP_to_netcdf/Tair/AWAP.Tair.3hr.2019.nc"
f = addfile (ncname,"r")

do day = 1,dom
if day .lt. 10 then
  fltname = "/g/data/w35/Shared_data/AWAP_raw_data/awap_tmax_C_day/bom-tmax-day-20090"+month+"0"+day+"-20090"+month+"0"+day+".flt"
else
  fltname = "/g/data/w35/Shared_data/AWAP_raw_data/awap_tmax_C_day/bom-tmax-day-20090"+month+day+"-20090"+month+day+".flt"
end if

tmp = fbindirread(fltname,0,(/681,841/),"float")
tmp = tmp(::-1,:)
tmp!0 ="lat"
tmp!1 ="lon"
tmp&lat = f->lat
tmp&lon = f->lon

res@tiMainString  = "Feb-"+day
plots = gsn_csm_contour_map(wks,tmp,res)

end do
end
