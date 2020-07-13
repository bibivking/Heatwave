;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

path = "/g/data/w35/mm3972/model/wrf/NUWRF/nuwrf_25km/LISWRF_configs/LIS_output"

var_name= (/\
            "Qle_tavg",\ ; "W m-2" "latent heat flux" ;
            "Qh_tavg",\  ; "W m-2" "sensible heat flux" ;
            "Qg_tavg",\  ; "W m-2" "soil heat flux" ;
            "Evap_tavg",\ ; "kg m-2 s-1" "total evapotranspiration" ;
            "TVeg_tavg",\ ; "kg m-2 s-1" "vegetation transpiration"
            "ESoil_tavg",\; "kg m-2 s-1" "bare soil evaporation"
            "Qs_tavg",\   ; "kg m-2 s-1" "surface runoff" ;
            "Qsb_tavg",\  ; "kg m-2 s-1" "subsurface runoff amount" ;
            "SnowDepth_inst",\; "m" "snow depth"
            "AvgSurfT_tavg",\; "K" "surface temperature" ;
            "SoilWet_inst",\ ; "-" "total soil wetness"
            "SoilTemp_inst",\;(time, 6, north_south, east_west) "K" "soil temperature" ;
            "SoilMoist_inst"\;(time, 6, north_south, east_west) "m^3 m-3" "soil moisture content"
           /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

;______________________________ Data _________________________________
year_s = 1979
year_e = 2008
year_sum = year_e-year_s+1

Qle      = new((/year_sum,139,159/),float)
Qh       = new((/year_sum,139,159/),float)
Qg       = new((/year_sum,139,159/),float)
Evap     = new((/year_sum,139,159/),float)
TVeg     = new((/year_sum,139,159/),float)
ESoil    = new((/year_sum,139,159/),float)
Qs       = new((/year_sum,139,159/),float)
Qsb      = new((/year_sum,139,159/),float)
SnowDepth= new((/year_sum,139,159/),float)
AvgSurfT = new((/year_sum,139,159/),float)
SoilWet  = new((/year_sum,139,159/),float)
SoilTemp = new((/year_sum,6,139,159/),float)
SoilMoist= new((/year_sum,6,139,159/),float)

Qle      = 0.0
Qh       = 0.0
Qg       = 0.0
Evap     = 0.0
TVeg     = 0.0
ESoil    = 0.0
Qs       = 0.0
Qsb      = 0.0
SnowDepth= 0.0
AvgSurfT = 0.0
SoilWet  = 0.0
SoilTemp = 0.0
SoilMoist= 0.0

do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    do month = 1, 12
      print("month = "+month)
      if month .lt. 10 then
        filename = path+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename = path+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      print(filename)
      f = addfile (filename,"r")

      Qle(year-year_s,:,:)        = Qle(year-year_s,:,:) + dim_avg_n_Wrap(f->Qle_tavg(:,:,:),0)/12
      Qh(year-year_s,:,:)         = Qh(year-year_s,:,:)  + dim_avg_n_Wrap(f->Qh_tavg(:,:,:),0)/12
      Qg(year-year_s,:,:)         = Qg(year-year_s,:,:)  + dim_avg_n_Wrap(f->Qg_tavg(:,:,:),0)/12
      Evap(year-year_s,:,:)       = Evap(year-year_s,:,:) + dim_avg_n_Wrap(f->Evap_tavg(:,:,:),0)*3600*24
      TVeg(year-year_s,:,:)       = TVeg(year-year_s,:,:) + dim_avg_n_Wrap(f->TVeg_tavg(:,:,:),0)*3600*24
      ESoil(year-year_s,:,:)      = ESoil(year-year_s,:,:)+ dim_avg_n_Wrap(f->ESoil_tavg(:,:,:),0)*3600*24
      Qs(year-year_s,:,:)         = Qs(year-year_s,:,:)   + dim_avg_n_Wrap(f->Qs_tavg(:,:,:),0)*3600*24
      Qsb(year-year_s,:,:)        = Qsb(year-year_s,:,:)  + dim_avg_n_Wrap(f->Qsb_tavg(:,:,:),0)*3600*24
      SnowDepth(year-year_s,:,:)  = SnowDepth(year-year_s,:,:) + dim_avg_n_Wrap(f->SnowDepth_inst(:,:,:),0)/12
      AvgSurfT(year-year_s,:,:)   = AvgSurfT(year-year_s,:,:)  + (dim_avg_n_Wrap(f->AvgSurfT_tavg(:,:,:),0)-273.15)/12
      SoilWet(year-year_s,:,:)    = SoilWet(year-year_s,:,:)   + dim_avg_n_Wrap(f->SoilWet_inst(:,:,:),0)/12
      SoilTemp(year-year_s,:,:,:) = SoilTemp(year-year_s,:,:,:)+ (dim_avg_n_Wrap(f->SoilTemp_inst(:,:,:,:),0)-273.15)/12
      SoilMoist(year-year_s,:,:,:)= SoilMoist(year-year_s,:,:,:)+ dim_avg_n_Wrap(f->SoilMoist_inst(:,:,:,:),0)/12
    end do ; month
end do ; year
print(Qle)
;_____________________________ PLOT ___________________________________
pic = "spitial_map"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

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
;
; A = -44.
; B = -10.
; C = 112.
; D = 154.
; res@trYMinF   = A
; res@trYMaxF   = B
; res@mpMinLatF = A
; res@mpMaxLatF = B
; res@trXMinF   = C
; res@trXMaxF   = D
; res@mpMinLonF = C
; res@mpMaxLonF = D
; res@gsnAddCyclic       = False
; res@tmYLValues=(/-10,-20,-30,-40/)
; res@tmYLLabels=(/"10~S~o~N~S","20~S~o~N~S","30~S~o~N~S","40~S~o~N~S"/)
; res@tmXBValues=(/120,130,140,150/)
; res@tmXBLabels=(/"120~S~o~N~E","130~S~o~N~E","140~S~o~N~E","150~S~o~N~E"/)

res@cnFillPalette = "WhiteBlueGreenYellowRed";"ViBlGrWhYeOrRe"
res@tiMainString = f->Qle_tavg@long_name+" ("+f->Qle_tavg@units+")"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(Qle,0),res)

res@tiMainString = f->Qh_tavg@long_name+" ("+f->Qh_tavg@units+")"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(Qh,0),res)

res@tiMainString = f->Qg_tavg@long_name+" ("+f->Qg_tavg@units+")"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(Qg,0),res)

res@tiMainString = f->Evap_tavg@long_name+" (mm/y)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(Evap,0),res)

res@tiMainString = f->TVeg_tavg@long_name+" (mm/y)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(TVeg,0),res)

res@tiMainString = f->ESoil_tavg@long_name+" (mm/y)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(ESoil,0),res)

res@tiMainString = f->Qs_tavg@long_name+" (mm/y)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(Qs,0),res)

res@tiMainString = f->Qsb_tavg@long_name+" (mm/y)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(Qsb,0),res)

res@tiMainString = f->SnowDepth_inst@long_name+" ("+f->SnowDepth_inst@units+")"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(SnowDepth,0),res)

res@tiMainString = f->AvgSurfT_tavg@long_name+" (C)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(AvgSurfT,0),res)

res@tiMainString = f->SoilWet_inst@long_name+" ("+f->SoilWet_inst@units+")"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(SoilWet,0),res)

res@tiMainString = f->SoilTemp_inst@long_name+" (C)"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(dim_avg_n_Wrap(SoilTemp,0),0),res)

res@tiMainString = f->SoilMoist_inst@long_name+" ("+f->SoilMoist_inst@units+")"
plot = gsn_csm_contour(wks,dim_avg_n_Wrap(dim_avg_n_Wrap(SoilMoist,0),0),res)

end
