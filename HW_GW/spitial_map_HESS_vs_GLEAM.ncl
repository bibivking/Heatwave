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

; =================================== Path =====================================
case_name   = (/"CTL-30x3+30yr/outputs-daily", "FREE_DRAIN/outputs-spinup30x3yr"/);, \
                ;"CTL-30x3+30yr-DeepRoot/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr-DeepRoot"/)
path_mum    = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_sum    = dimsizes(case_name)
path        = new(case_sum,"string")
pic_message = "CTL_vs_FD"
order_name  = (/"(c)","(d)","(e)","(f)","(g)","(h)"/) ; "(a)","(b)",

do case_num = 0, case_sum -1
   path(case_num) = path_mum+case_name(case_num);+"/outputs"
end do

path_GLEAM= "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_3/3_3a/yearly";

; ================================= Predefine ==================================
year_s   = 2000
year_e   = 2019
year_sum = year_e-year_s+1
row_sum  = 2
col_sum  = 3

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; ================================== OPTIONS ===================================
plot_type = "all"
scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"

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

filename  = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f         = addfile(filename,"r")

lat_sum   = dimsizes(f->Rainf(0,{A:B},0))
lon_sum   = dimsizes(f->Rainf(0,0,{C:D}))
lat_size  = dimsizes(f->y)
lon_size  = dimsizes(f->x)
CABLE_lat = f->y
CABLE_lon = f->x
default_val = -1.e+33
delete(f)

Evap      = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val)
; Fwsoil    = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val)
; Qh        = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val)
; VegT      = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val)
; Tair      = new((/case_sum, year_sum, lat_sum, lon_sum/),float, default_val)

Evap      = 0.
; Fwsoil    = 0.
; Qh        = 0.
; VegT      = 0.
; Tair      = 0.

do case_num = 0, case_sum -1
   do year = year_s, year_e

      print("year = " +year)

      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if

      filename = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
      f        = addfile (filename,"r")

      Evap(case_num,year-year_s,:,:)     = dim_sum_n_Wrap(f->Evap(:,{A:B},{C:D}),0)
      ; Fwsoil(case_num,year-year_s,:,:)   = dim_avg_n_Wrap(f->Fwsoil(:,{A:B},{C:D}),0)
      ; Qh(case_num,year-year_s,:,:)       = dim_avg_n_Wrap(f->Qh(:,{A:B},{C:D}),0)
      ; VegT(case_num,year-year_s,:,:)     = dim_avg_n_Wrap(f->VegT(:,{A:B},{C:D}),0)
      ; Tair(case_num,year-year_s,:,:)     = dim_avg_n_Wrap(f->Tair(:,{A:B},{C:D}),0)
      delete(f)
   end do ; year
end do ; case_num

Evap     = Evap*3600.*24.

; VegT     = VegT-273.15
; Tair     = Tair-273.15

; ================================ READ GLEAM ==================================
; use 2000-2019
gleam_ET = path_GLEAM+"/E_1980_2018_GLEAM_v3.3a_YR.nc"
G_ET     = addfile (gleam_ET,"r")

E_tmp    = G_ET->E(:,{C:D},{A:B})

Evap_GLEAM = E_tmp(time|20:38,lat|:, lon|:)
printVarSummary(Evap_GLEAM)
delete(G_ET)

; ___________________________ Regrid GLEAM ______________________________
Evap_GLEAM_regrid = linint2_Wrap(Evap_GLEAM&lon,Evap_GLEAM&lat,Evap_GLEAM,False,CABLE_lon,CABLE_lat,0)
printVarSummary(Evap_GLEAM_regrid)


; ================================= Plotting ===================================
res                         = True
res@gsnDraw                 = False  ; Don't draw plots
res@gsnFrame                = False  ; ����ҳ
res@gsnMaximize             = True

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

res@lbLabelBarOn           = True
res@lbBoxLinesOn           = True                       ;�ر�lbar box ����
res@lbTitleFont            = 25
res@lbLabelFont            = 25
res@lbTitleFontHeightF     = 0.013
res@lbLabelFontHeightF     = 0.013
; res@lbOrientation           = "Vertical"
res@lbTitlePosition        = "Bottom"
res@pmLabelBarHeightF       = 0.1
res@pmLabelBarParallelPosF  = 0.5
res@pmLabelBarOrthogonalPosF= 0.09

res@tmXBLabelFont           = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont           = 25
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF    = 0.015
res@tmYLLabelFontHeightF    = 0.015
; res@mpDataBaseVersion       = "HighRes"

res@tmYLMode                = "Explicit"
res@tmXBMode                = "Explicit"
res@txString                = ""
res@gsnLeftString           = ""
res@gsnRightString          = ""

if scale .eq. "Global" then
  res@trYMinF         = -60
  res@trYMaxF         = 80
  res@mpMinLatF       = -60
  res@mpMaxLatF       = 80
  res@mpCenterLonF    = 160
  res@tmYLMode        = "Explicit"
  res@tmYLValues      =(/60,30,0,-30,-60/)
  res@tmYLLabels      =(/"60~S~o~N~N","30~S~o~N~N","0~S~o~N~","30~S~o~N~S","60~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      =(/240,300,0,60,120,180/)
  res@tmXBLabels      =(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)
else if scale .eq. "SE-AU" then
  res@gsnAddCyclic    = False
  res@trYMinF         = -40.
  res@trYMaxF         = -28.
  res@mpMinLatF       = -40.
  res@mpMaxLatF       = -28.
  res@trXMinF         = 140.
  res@trXMaxF         = 154.
  res@mpMinLonF       = 140.
  res@mpMaxLonF       = 154.
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-40,-35,-30,-25/)
  res@tmYLLabels      = (/"40~S~o~N~S","35~S~o~N~S","30~S~o~N~S","25~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/140,145,150/)
  res@tmXBLabels      = (/"140~S~o~N~E","145~S~o~N~E","150~S~o~N~E"/)
else if scale .eq. "AU" then
  res@gsnAddCyclic    = False
  res@trYMinF         = -44
  res@trYMaxF         = -10
  res@mpMinLatF       = -44
  res@mpMaxLatF       = -10
  res@trXMinF         = 112
  res@trXMaxF         = 154
  res@mpMinLonF       = 112
  res@mpMaxLonF       = 154
  res@tmYLMode        = "Explicit"
  res@tmYLValues      = (/-40,-35,-30,-25,-20,-15,-10/)
  res@tmYLLabels      = (/"40~S~o~N~N","35~S~o~N~N","30~S~o~N~","25~S~o~N~S","20~S~o~N~N","15~S~o~N~N","10~S~o~N~"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/115,120,125,130,135,140,145,150/)
  res@tmXBLabels      = (/"115~S~o~N~W","120~S~o~N~W","125~S~o~N~","130~S~o~N~W","135~S~o~N~W","140~S~o~N~","145~S~o~N~W","150~S~o~N~"/)
end if
end if
end if

;**************** ͼ�������� *****************
pres                    = True
pres@wkWidth            = 2500
pres@wkHeight           = 2500

pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False
pres@cnLevelSelectionMode = "ExplicitLevels"
pres@lbBoxLinesOn       = True
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01

txres                   = True
txres@txPerimOn         = False
txres@txFont            = 25
txres@txFontHeightF     = 0.05

amres                   = True
amres@amParallelPosF    = -0.485    ; This is the right edge of the plot.
amres@amOrthogonalPosF  = -0.46    ; This is the bottom edge of the plot.
amres@amJust            = "TopLeft"

; =========================== plot Evaptranspiration ===========================
plots = new((case_sum+1)*2,graphic)

pic = "./plots/HESS_vs_GLEAM_"+pic_message
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

counter = 0

; ___ 2001-2009 ___
res@lbLabelBarOn = True

; GLEAM -----
res@cnLevels  = fspan(0,1000,21)
res@cnFillPalette = "vegetation_ClarkU";"GMT_seis";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
res@lbTitleString = "E (mm yr~S~-1~N~)"

plots(counter) = gsn_csm_contour_map(wks,dim_avg_n_Wrap(Evap_GLEAM(1:9,:,:),0),res)
txid1   = gsn_create_text(wks, order_name(counter), txres)
annoid1 = gsn_add_annotation(plots(counter), txid1, amres)
counter = counter+1
delete(res@cnLevels)

; -----
res@cnLevels  = (/-400,-350,-300,-250,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400/)
delete(res@cnFillPalette)
cmap = read_colormap_file("BlAqGrWh2YeOrReVi22");("NEO_div_vegetation_a")("ViBlGrWhYeOrRe");
res@cnFillPalette = cmap(::-1,:)
res@lbTitleString = "~F5~D"+"~F25~E (mm yr~S~-1~N~)"

; CTL - GLEAM
tmp_GLEAM = doubletofloat(dim_avg_n_Wrap(Evap_GLEAM_regrid(1:9,:,:),0))
tmp       = dim_avg_n_Wrap(Evap(0,1:9,:,:),0)
tmp       = tmp - tmp_GLEAM

printVarSummary(tmp)
plots(counter) = gsn_csm_contour_map(wks, tmp, res)
txid1          = gsn_create_text(wks, order_name(counter), txres)
annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
counter        = counter + 1
delete(tmp)

; FD - GLEAM
tmp = dim_avg_n_Wrap(Evap(1,1:9,:,:),0)
tmp = tmp - tmp_GLEAM
plots(counter) = gsn_csm_contour_map(wks, tmp, res)
txid1          = gsn_create_text(wks, order_name(counter), txres)
annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
counter        = counter + 1
delete(res@cnLevels)
delete(tmp)
delete(tmp_GLEAM)
;
; ; Fwsoil FREE_DRAIN - CTL
; res@cnLevels  = fspan(-0.4,0.0,21)
; res@cnFillPalette = "MPL_hot";"MPL_autumn";"MPL_BrBG";"ViBlGrWhYeOrRe"
;
; tmp = Fwsoil(1,1:9,:,:)
; tmp = tmp - Fwsoil(0,1:9,:,:)
; printVarSummary(tmp)
; plots(counter)  = gsn_csm_contour_map(wks, dim_avg_n_Wrap(tmp,0),res)
; txid1          = gsn_create_text(wks, order_name(counter), txres)
; annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
; counter        = counter + 1
; delete(res@cnLevels)
; delete(tmp)
;
; ; ___ 2017-2018 ___
res@lbLabelBarOn            = True

; GLEAM 2017-2018
delete(res@cnFillPalette)
res@cnLevels  = fspan(0,1000,21)
res@cnFillPalette = "vegetation_ClarkU";"GMT_seis";"NOC_ndvi";"MPL_RdYlGn";"WhiteBlueGreenYellowRed"
res@lbTitleString = "E (mm yr~S~-1~N~)"

; tmp = dim_avg_n_Wrap(Evap_GLEAM(17:18,:,:),0)
; tmp = tmp - dim_avg_n_Wrap(Evap_GLEAM(1:9,:,:),0)
plots(counter) = gsn_csm_contour_map(wks,dim_avg_n_Wrap(Evap_GLEAM(17:18,:,:),0),res)
txid1   = gsn_create_text(wks, order_name(counter), txres)
annoid1 = gsn_add_annotation(plots(counter), txid1, amres)
counter = counter+1
delete(res@cnLevels)
; delete(tmp)

; ----
tmp_GLEAM     = doubletofloat(dim_avg_n_Wrap(Evap_GLEAM_regrid(17:18,:,:),0))
res@cnLevels  = (/-400,-350,-300,-250,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400/)

delete(res@cnFillPalette)
cmap = read_colormap_file("BlAqGrWh2YeOrReVi22");("NEO_div_vegetation_a")("ViBlGrWhYeOrRe");
res@cnFillPalette = cmap(::-1,:)

res@lbTitleString = "~F5~D"+"~F25~E (mm yr~S~-1~N~)"

; CTL - GLEAM
tmp = dim_avg_n_Wrap(Evap(0,17:18,:,:),0)
tmp = tmp - tmp_GLEAM ; dim_avg_n_Wrap(Evap(0,1:9,:,:),0) ;
printVarSummary(tmp)
plots(counter) = gsn_csm_contour_map(wks, tmp, res)
txid1          = gsn_create_text(wks, order_name(counter), txres)
annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
counter        = counter + 1
delete(tmp)

; FD - GLEAM
tmp = dim_avg_n_Wrap(Evap(1,17:18,:,:),0)
tmp = tmp - tmp_GLEAM ;dim_avg_n_Wrap(Evap(1,1:9,:,:),0)
plots(counter) = gsn_csm_contour_map(wks, tmp, res)
txid1          = gsn_create_text(wks, order_name(counter), txres)
annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
counter        = counter + 1
delete(res@cnLevels)
delete(tmp)
; delete(tmp_GLEAM)

gsn_panel(wks,(/plots/),(/row_sum,col_sum/),pres)

end

; ; =========================== plot deltaT & Qh ===========================
; plots = new(case_sum*4,graphic)
;
; pic = "./plots/HESS_vs_GLEAM_Qh_"+pic_message
; wks = gsn_open_wks("pdf",pic)
; gsn_define_colormap(wks,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")
;
; counter = 0
;
; ; ___ 2001-2009 ___
; ; delta T CTL
; res@cnLevels  = fspan(0,4,21)
; res@cnFillPalette = "WhiteBlueGreenYellowRed"
;
; tmp = VegT(0,1:9,:,:)
; tmp = tmp - Tair(0,1:9,:,:)
; plots(counter) = gsn_csm_contour_map(wks,dim_avg_n_Wrap(tmp,0),res)
; txid1   = gsn_create_text(wks, order_name(counter), txres)
; annoid1 = gsn_add_annotation(plots(counter), txid1, amres)
; counter = counter+1
; delete(tmp)
;
; ; delta T FREE_DRAIN
; res@cnLevels  = fspan(-1,1,21)
; res@cnFillPalette = "ViBlGrWhYeOrRe"
;
; tmp = VegT(1,1:9,:,:)
; tmp = tmp - Tair(1,1:9,:,:) - (VegT(0,1:9,:,:) - Tair(0,1:9,:,:))
; plots(counter)  = gsn_csm_contour_map(wks, dim_avg_n_Wrap(tmp,0),res)
; txid1          = gsn_create_text(wks, order_name(counter), txres)
; annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
; counter        = counter + 1
; delete(res@cnLevels)
; delete(tmp)
;
; ; Qh CTL
; res@cnLevels  = fspan(0,160,17)
; res@cnFillPalette = "WhiteBlueGreenYellowRed"
; plots(counter)  = gsn_csm_contour_map(wks, dim_avg_n_Wrap(Qh(0,1:9,:,:),0),res)
; txid1          = gsn_create_text(wks, order_name(counter), txres)
; annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
; counter        = counter + 1
; delete(res@cnLevels)
;
; ; Qh FREE_DRAIN
; res@cnLevels  = fspan(-20,20,11)
; res@cnFillPalette = "ViBlGrWhYeOrRe"
; tmp = Qh(1,1:9,:,:)
; tmp = tmp - Qh(0,1:9,:,:)
; plots(counter)  = gsn_csm_contour_map(wks, dim_avg_n_Wrap(tmp,0),res)
; txid1          = gsn_create_text(wks, order_name(counter), txres)
; annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
; counter        = counter + 1
; delete(res@cnLevels)
; delete(tmp)
;
; ; ___ 2017-2019 ___
; ; delta T CTL
; res@cnLevels  = fspan(0,4,21)
; res@cnFillPalette = "WhiteBlueGreenYellowRed"
;
; tmp = VegT(0,16:18,:,:)
; tmp = tmp - Tair(0,16:18,:,:)
; plots(counter) = gsn_csm_contour_map(wks,dim_avg_n_Wrap(tmp,0),res)
; txid1   = gsn_create_text(wks, order_name(counter), txres)
; annoid1 = gsn_add_annotation(plots(counter), txid1, amres)
; counter = counter+1
; delete(tmp)
;
; ; delta T FREE_DRAIN
; res@cnLevels  = fspan(-1,1,21)
; res@cnFillPalette = "ViBlGrWhYeOrRe"
;
; tmp = VegT(1,16:18,:,:)
; tmp = tmp - Tair(1,16:18,:,:) - (VegT(0,16:18,:,:) - Tair(0,16:18,:,:))
; plots(counter)  = gsn_csm_contour_map(wks, dim_avg_n_Wrap(tmp,0),res)
; txid1          = gsn_create_text(wks, order_name(counter), txres)
; annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
; counter        = counter + 1
; delete(res@cnLevels)
; delete(tmp)
;
; ; Qh CTL
; res@cnLevels  = fspan(0,160,17)
; res@cnFillPalette = "WhiteBlueGreenYellowRed"
; plots(counter)  = gsn_csm_contour_map(wks, dim_avg_n_Wrap(Qh(0,16:18,:,:),0),res)
; txid1          = gsn_create_text(wks, order_name(counter), txres)
; annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
; counter        = counter + 1
; delete(res@cnLevels)
;
; ; Qh FREE_DRAIN
; res@cnLevels  = fspan(-20,20,11)
; res@cnFillPalette = "ViBlGrWhYeOrRe"
; tmp = Qh(1,16:18,:,:)
; tmp = tmp - Qh(0,16:18,:,:)
; plots(counter)  = gsn_csm_contour_map(wks, dim_avg_n_Wrap(tmp,0),res)
; txid1          = gsn_create_text(wks, order_name(counter), txres)
; annoid1        = gsn_add_annotation(plots(counter), txid1, amres)
; counter        = counter + 1
; delete(res@cnLevels)
; delete(tmp)
