;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 500000000000 ;
end setvalues
;_____________________________________


begin

path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ERAI_check_hydsx10/LIS_output_GWwb03"
path2 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ERAI_check_hydsx10/LIS_output_GWwb045"
path3 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ERAI_check_GWextrct/LIS_output_hydsx1"
path4 = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/ERAI_check_GWextrct/LIS_output_hydsx10"

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; ____________________________ OPTIONS _______________________________
hydr_param = True
scale      = "domain"
case_name  = "GWwb_03_vs_045_vs_hyds10"

; region: AU
A = -44.
B = -10.
C = 112.
D = 154.

ilat = 15
ilon = 30

;______________________________ Data _________________________________
year_s = 2018
year_e = 2018
year_sum = year_e-year_s+1
time_tot = year_sum*365*48 ;365

Rain     = new((/4,time_tot/),float)
Wind     = new((/4,time_tot/),float)
Qair     = new((/4,time_tot/),float)
Tair     = new((/4,time_tot/),float)
Psurf    = new((/4,time_tot/),float)
SWdown   = new((/4,time_tot/),float)
LWdown   = new((/4,time_tot/),float)

Qle      = new((/4,time_tot/),float)
Qh       = new((/4,time_tot/),float)
Qg       = new((/4,time_tot/),float)
Evap     = new((/4,time_tot/),float)
TVeg     = new((/4,time_tot/),float)
ESoil    = new((/4,time_tot/),float)
Qs       = new((/4,time_tot/),float)
Qsb      = new((/4,time_tot/),float)
Qrec     = new((/4,time_tot/),float)
; AvgSurfT = new((/4,time_tot/),float)

if hydr_param then
  GWwb     = new((/4,time_tot/),float)
  ; Fwsoil   = new((/4,time_tot/),float)
  WatTab   = new((/4,time_tot/),float)
else
  SnowDepth= new((/4,time_tot/),float)
end if

; SoilTemp = new((/4,time_tot,6/),float)
SoilMoist= new((/4,time_tot,6/),float)

Rain     = 0.
Wind     = 0.
Qair     = 0.
Tair     = 0.
Psurf    = 0.
SWdown   = 0.
LWdown   = 0.
Qle      = 0.
Qh       = 0.
Qg       = 0.
Evap     = 0.
TVeg     = 0.
ESoil    = 0.
Qs       = 0.
Qsb      = 0.
Qrec     = 0.
; AvgSurfT = 0.

if hydr_param then
  GWwb     = 0.
  ; Fwsoil   = 0.
  WatTab   = 0.
else
  SnowDepth= 0.
end if

; SoilTemp = 0.
SoilMoist= 0.


; ------ print paramteres -------

filename = path1+"/LIS.CABLE.2018090100.d01.nc"
f = addfile (filename,"r")

print("lat : "+ f->lat(1,ilat,ilon) )
print("lon : "+ f->lon(1,ilat,ilon) )
print("soil type : "+ f->Soiltype_inst(1,ilat,ilon) )
print("veg type : "+ f->Landcover_inst(1,ilat,ilon) )
; print("bch : "+ f->bch_inst(1,ilat,ilon) )
; print("sfc : "+ f->sfc_inst(1,ilat,ilon) )
; print("swilt : "+ f->swilt_inst(1,ilat,ilon) )
; print("ssat : "+ f->Porosity_inst(1,ilat,ilon) )
; print("hyds : "+ f->hyds_inst(1,ilat,ilon) )
; print("sucs : "+ f->sucs_inst(1,ilat,ilon) )
SWdown_map = f->SWdown_f_inst(1,:,:)
SWdown_map = 0.

i = 0
do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       tsom = (/1488, 1392,1488,1440,1488,1440,1488,1488,1440,1488,1440,1488/)
    else
       tsom = (/1488, 1344,1488,1440,1488,1440,1488,1488,1440,1488,1440,1488/)
    end if

    ; read LIS-CABLE
    do month = 1, 12
      print("month = "+month)
      if month .lt. 10 then
        filename = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
        filename1 = path2+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
        filename2 = path3+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
        filename3 = path4+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
        filename1 = path2+"/LIS.CABLE."+year+month+"0100.d01.nc"
        filename2 = path3+"/LIS.CABLE."+year+month+"0100.d01.nc"
        filename3 = path4+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      f  = addfile (filename,"r")
      f1 = addfile (filename1,"r")
      f2 = addfile (filename2,"r")
      f3 = addfile (filename3,"r")

      ; ========================= READING =============================

      do ts = 1, tsom(month-1)

      ; =================== Met Force ====================
        Rain(0,i)     = avg(f->Rainf_f_inst(ts-1,ilat,ilon))*3600*0.5
        Wind(0,i)     = avg(f->Wind_f_inst(ts-1,ilat,ilon))
        Qair(0,i)     = avg(f->Qair_f_inst(ts-1,ilat,ilon))
        Tair(0,i)     = avg(f->Tair_f_inst(ts-1,ilat,ilon))-273.15
        Psurf(0,i)    = avg(f->Psurf_f_inst(ts-1,ilat,ilon))
        SWdown(0,i)   = avg(f->SWdown_f_inst(ts-1,ilat,ilon))
        LWdown(0,i)   = avg(f->LWdown_f_inst(ts-1,ilat,ilon))
        SWdown_map    = SWdown_map + f->SWdown_f_inst(ts-1,:,:)


        Rain(1,i)     = avg(f1->Rainf_f_inst(ts-1,ilat,ilon))*3600*0.5
        Wind(1,i)     = avg(f1->Wind_f_inst(ts-1,ilat,ilon))
        Qair(1,i)     = avg(f1->Qair_f_inst(ts-1,ilat,ilon))
        Tair(1,i)     = avg(f1->Tair_f_inst(ts-1,ilat,ilon))-273.15
        Psurf(1,i)    = avg(f1->Psurf_f_inst(ts-1,ilat,ilon))
        SWdown(1,i)   = avg(f1->SWdown_f_inst(ts-1,ilat,ilon))
        LWdown(1,i)   = avg(f1->LWdown_f_inst(ts-1,ilat,ilon))

        Rain(2,i)     = avg(f2->Rainf_f_inst(ts-1,ilat,ilon))*3600*0.5
        Wind(2,i)     = avg(f2->Wind_f_inst(ts-1,ilat,ilon))
        Qair(2,i)     = avg(f2->Qair_f_inst(ts-1,ilat,ilon))
        Tair(2,i)     = avg(f2->Tair_f_inst(ts-1,ilat,ilon))-273.15
        Psurf(2,i)    = avg(f2->Psurf_f_inst(ts-1,ilat,ilon))
        SWdown(2,i)   = avg(f2->SWdown_f_inst(ts-1,ilat,ilon))
        LWdown(2,i)   = avg(f2->LWdown_f_inst(ts-1,ilat,ilon))


        Rain(3,i)     = avg(f3->Rainf_f_inst(ts-1,ilat,ilon))*3600*0.5
        Wind(3,i)     = avg(f3->Wind_f_inst(ts-1,ilat,ilon))
        Qair(3,i)     = avg(f3->Qair_f_inst(ts-1,ilat,ilon))
        Tair(3,i)     = avg(f3->Tair_f_inst(ts-1,ilat,ilon))-273.15
        Psurf(3,i)    = avg(f3->Psurf_f_inst(ts-1,ilat,ilon))
        SWdown(3,i)   = avg(f3->SWdown_f_inst(ts-1,ilat,ilon))
        LWdown(3,i)   = avg(f3->LWdown_f_inst(ts-1,ilat,ilon))

      ; =================== Variables ===================
        Qle(0,i)      = avg(f->Qle_tavg(ts-1,ilat,ilon))
        Qh(0,i)       = avg(f->Qh_tavg(ts-1,ilat,ilon))
        Qg(0,i)       = avg(f->Qg_tavg(ts-1,ilat,ilon))
        Evap(0,i)     = avg(f->Evap_tavg(ts-1,ilat,ilon))*3600*0.5
        TVeg(0,i)     = avg(f->TVeg_tavg(ts-1,ilat,ilon))*3600*0.5
        ESoil(0,i)    = avg(f->ESoil_tavg(ts-1,ilat,ilon))*3600*0.5
        Qs(0,i)       = avg(f->Qs_tavg(ts-1,ilat,ilon))*3600*0.5
        Qsb(0,i)      = avg(f->Qsb_tavg(ts-1,ilat,ilon))*3600*0.5
        Qrec(0,i)     = avg(f->Qrec_tavg(ts-1,ilat,ilon))*3600*0.5

        Qle(1,i)      = avg(f1->Qle_tavg(ts-1,ilat,ilon))
        Qh(1,i)       = avg(f1->Qh_tavg(ts-1,ilat,ilon))
        Qg(1,i)       = avg(f1->Qg_tavg(ts-1,ilat,ilon))
        Evap(1,i)     = avg(f1->Evap_tavg(ts-1,ilat,ilon))*3600*0.5
        TVeg(1,i)     = avg(f1->TVeg_tavg(ts-1,ilat,ilon))*3600*0.5
        ESoil(1,i)    = avg(f1->ESoil_tavg(ts-1,ilat,ilon))*3600*0.5
        Qs(1,i)       = avg(f1->Qs_tavg(ts-1,ilat,ilon))*3600*0.5
        Qsb(1,i)      = avg(f1->Qsb_tavg(ts-1,ilat,ilon))*3600*0.5
        Qrec(1,i)     = avg(f1->Qrec_tavg(ts-1,ilat,ilon))*3600*0.5

        Qle(2,i)      = avg(f2->Qle_tavg(ts-1,ilat,ilon))
        Qh(2,i)       = avg(f2->Qh_tavg(ts-1,ilat,ilon))
        Qg(2,i)       = avg(f2->Qg_tavg(ts-1,ilat,ilon))
        Evap(2,i)     = avg(f2->Evap_tavg(ts-1,ilat,ilon))*3600*0.5
        TVeg(2,i)     = avg(f2->TVeg_tavg(ts-1,ilat,ilon))*3600*0.5
        ESoil(2,i)    = avg(f2->ESoil_tavg(ts-1,ilat,ilon))*3600*0.5
        Qs(2,i)       = avg(f2->Qs_tavg(ts-1,ilat,ilon))*3600*0.5
        Qsb(2,i)      = avg(f2->Qsb_tavg(ts-1,ilat,ilon))*3600*0.5
        Qrec(2,i)     = avg(f2->Qrec_tavg(ts-1,ilat,ilon))*3600*0.5

        Qle(3,i)      = avg(f3->Qle_tavg(ts-1,ilat,ilon))
        Qh(3,i)       = avg(f3->Qh_tavg(ts-1,ilat,ilon))
        Qg(3,i)       = avg(f3->Qg_tavg(ts-1,ilat,ilon))
        Evap(3,i)     = avg(f3->Evap_tavg(ts-1,ilat,ilon))*3600*0.5
        TVeg(3,i)     = avg(f3->TVeg_tavg(ts-1,ilat,ilon))*3600*0.5
        ESoil(3,i)    = avg(f3->ESoil_tavg(ts-1,ilat,ilon))*3600*0.5
        Qs(3,i)       = avg(f3->Qs_tavg(ts-1,ilat,ilon))*3600*0.5
        Qsb(3,i)      = avg(f3->Qsb_tavg(ts-1,ilat,ilon))*3600*0.5
        Qrec(3,i)     = avg(f3->Qrec_tavg(ts-1,ilat,ilon))*3600*0.5

        ; AvgSurfT(0,i) = avg(f->AvgSurfT_tavg(ts-1,ilat,ilon))-273.15
        ; AvgSurfT(1,i) = avg(f1->AvgSurfT_tavg(ts-1,ilat,ilon))-273.15

        if hydr_param then
          GWwb(0,i)    = avg(f->GWwb_inst(ts-1,ilat,ilon))
          ; Fwsoil(0,i)  = avg(f->fwsoil_inst(ts-1,ilat,ilon))
          WatTab(0,i)  = avg(f->WaterTableD_tavg(ts-1,ilat,ilon))/1000.


          GWwb(1,i)    = avg(f1->GWwb_inst(ts-1,ilat,ilon))
          ; Fwsoil(1,i)  = avg(f1->fwsoil_inst(ts-1,ilat,ilon))
          WatTab(1,i)  = avg(f1->WaterTableD_tavg(ts-1,ilat,ilon))/1000.

          GWwb(2,i)    = avg(f2->GWwb_inst(ts-1,ilat,ilon))
          ; Fwsoil(2,i)  = avg(f2->fwsoil_inst(ts-1,ilat,ilon))
          WatTab(2,i)  = avg(f2->WaterTableD_tavg(ts-1,ilat,ilon))/1000.

          GWwb(3,i)    = avg(f3->GWwb_inst(ts-1,ilat,ilon))
          ; Fwsoil(3,i)  = avg(f3->fwsoil_inst(ts-1,ilat,ilon))
          WatTab(3,i)  = avg(f3->WaterTableD_tavg(ts-1,ilat,ilon))/1000.
        else
          SnowDepth(0,i)= avg(f->SnowDepth_inst(ts-1,ilat,ilon))
          SnowDepth(1,i)= avg(f1->SnowDepth_inst(ts-1,ilat,ilon))
        end if

        do soil = 0,5
          ; SoilTemp(0,i,soil) = avg(f->SoilTemp_inst(ts-1,soil,ilat,ilon))-273.15
          SoilMoist(0,i,soil)= avg(f->SoilMoist_inst(ts-1,soil,ilat,ilon))
          ; SoilTemp(1,i,soil) = avg(f1->SoilTemp_inst(ts-1,soil,ilat,ilon))-273.15
          SoilMoist(1,i,soil)= avg(f1->SoilMoist_inst(ts-1,soil,ilat,ilon))
          ; SoilTemp(2,i,soil) = avg(f2->SoilTemp_inst(ts-1,soil,ilat,ilon))-273.15
          SoilMoist(2,i,soil)= avg(f2->SoilMoist_inst(ts-1,soil,ilat,ilon))
          ; SoilTemp(3,i,soil) = avg(f3->SoilTemp_inst(ts-1,soil,ilat,ilon))-273.15
          SoilMoist(3,i,soil)= avg(f3->SoilMoist_inst(ts-1,soil,ilat,ilon))
        end do
        i = i + 1
      end do ; ts
      ; delete(filename)
      ; delete(f)
    end do ; month
    ;delete(f1)
end do ; year
print("prepare to plot")
print(avg(SWdown(0,:)))

;_____________________________ PLOT ___________________________________
if hydr_param then
  pic = "./plots/Spinup_LIS_"+case_name+"_"+scale+"_"+year_s+"-"+year_e
else
  pic = "./plots/Spinup_LIS_"+case_name+"_"+scale+"_"+year_s+"-"+year_e
end if

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��
res                    = True
res@gsnMaximize        = True
res@vpHeightF          = 0.6
res@vpWidthF           = 0.6
res@vpXF               = 0.3
res@vpYF               = 0.85
res@pmLegendWidthF     = 0.14                     ; Resize legend width
res@pmLegendHeightF    = 0.11                     ; and height

res@xyLineThicknesses  = (/  1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0/)          ; make second line thicker
res@pmLegendDisplayMode = "Always"
res@xyLineColors       = (/"blue","red","green","orange","red","red","red"/)          ; change line color
res@xyExplicitLegendLabels = (/"GWwb=0.3 hydsx10","GWwb=0.45 hydsx10","GWextrct hydsx1","GWextrct hydsx10"/)
;res@trYMinF = 0


x             = ispan(1,time_tot,1)
res@tmXBMode  = "Explicit"
res@tmXBValues= ispan(1,time_tot+1,10*48)
res@tmXBLabels= tostring(ispan(1,365+1,10)) ; 365 ;365


filename = path1+"/LIS.CABLE.2018090100.d01.nc"
f = addfile(filename,"r")


; SWdown

res1                    = True
res1@cnFillMode         = "RasterFill"            ; Raster Mode
res1@cnFillOn           = True                            ;��ɫ
res1@tmBorderThicknessF = 3.0
res1@gsnDraw            = True  ; Don't draw plots
res1@gsnFrame           = True  ; ����ҳ
res1@lbLabelBarOn       = True
res1@lbBoxLinesOn       = True
res1@cnLinesOn          = False
res1@cnLineLabelsOn     = False
res1@gsnMaximize        = True
res1@cnExplicitLabelBarLabelsOn = True   ;?
res1@cnLevelSelectionMode = "ExplicitLevels"
res1@tiMainString = "SWdown (W/m^2)"
levels = fspan(0,500,51)
res1@cnLevels = levels
plot = gsn_csm_contour(wks,SWdown_map/time_tot,res1)
delete(levels)
delete(res1@cnLevels)


print("plot rainfall")
; Rain_rc    = regline(x, Rain)          ; slope
; Rain_arry  = new ( (/2,time_tot/), typeof(Rain))
; Rain_arry(0,:) = Rain
; Rain_arry(1,:) = Rain_rc*x + Rain_rc@yintercept
if hydr_param then
  res@tiYAxisString  = f->TotalPrecip_tavg@long_name+" (mm/0.5hr)"
else
  res@tiYAxisString  = f->Rainf_f_inst@long_name+" (mm/0.5hr)"
end if
res@tiXAxisString  = "tss"
plot = gsn_csm_xy(wks,x,Rain,res)

if hydr_param then
  res@tiYAxisString  = "Wind Speed (m/s)"
else
  res@tiYAxisString  = "Wind Speed (m/s)"
end if
res@tiXAxisString  = "tss"
plot = gsn_csm_xy(wks,x,Wind,res)

if hydr_param then
  res@tiYAxisString  = "Specific Humidity (kg/kg)"
else
  res@tiYAxisString  = "Specific Humidity (kg/kg)"
end if
res@tiXAxisString  = "tss"
plot = gsn_csm_xy(wks,x,Qair,res)

if hydr_param then
  res@tiYAxisString  = "Air temperature (degree C)"
else
  res@tiYAxisString  = "Air temperature (degree C)"
end if
res@tiXAxisString  = "tss"
plot = gsn_csm_xy(wks,x,Tair,res)

if hydr_param then
  res@tiYAxisString  = "Air pressure (Pa)"
else
  res@tiYAxisString  = "Air pressure (Pa)"
end if
res@tiXAxisString  = "tss"
plot = gsn_csm_xy(wks,x,Psurf,res)

if hydr_param then
  res@tiYAxisString  = "Downwelling Shortwave Flux (W m-2)"
else
  res@tiYAxisString  = "Downwelling Shortwave Flux(W m-2)"
end if
res@tiXAxisString  = "tss"
plot = gsn_csm_xy(wks,x,SWdown,res)

if hydr_param then
  res@tiYAxisString  = "Downwelling Longwave Flux (W m-2)"
else
  res@tiYAxisString  = "Downwelling Longwave Flux (W m-2)"
end if
res@tiXAxisString  = "tss"
plot = gsn_csm_xy(wks,x,LWdown,res)

print("plot latent heat")
; Qle_rc    = regline(x, Qle)          ; slope
; Qle_arry  = new ( (/2,time_tot/), typeof(Qle))
; Qle_arry(0,:) = Qle
; Qle_arry(1,:) = Qle_rc*x + Qle_rc@yintercept
res@tiYAxisString  = f->Qle_tavg@long_name+" ("+f->Qle_tavg@units+")"
res@tiXAxisString  = "tss"
plot = gsn_csm_xy(wks,x,Qle,res)

print("plot sensible heat")
; Qh_rc    = regline(x, Qh)          ; slope
; Qh_arry  = new ( (/2,time_tot/), typeof(Qh))
; Qh_arry(0,:) = Qh
; Qh_arry(1,:) = Qh_rc*x + Qh_rc@yintercept
res@tiYAxisString  = f->Qh_tavg@long_name+" ("+f->Qh_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qh,res)

print("plot ground heat")
; Qg_rc    = regline(x, Qg)          ; slope
; Qg_arry  = new ( (/2,time_tot/), typeof(Qg))
; Qg_arry(0,:) = Qg
; Qg_arry(1,:) = Qg_rc*x + Qg_rc@yintercept
res@tiYAxisString  = f->Qg_tavg@long_name+" ("+f->Qg_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qg,res)

print("plot ET")
; Evap_rc    = regline(x, Evap)          ; slope
; Evap_arry  = new ( (/2,time_tot/), typeof(Evap))
; Evap_arry(0,:) = Evap
; Evap_arry(1,:) = Evap_rc*x + Evap_rc@yintercept
res@tiYAxisString  = f->Evap_tavg@long_name+" (mm/0.5hr)"
plot = gsn_csm_xy(wks,x,Evap,res)

print("plot trans")
; TVeg_rc    = regline(x, TVeg)          ; slope
; TVeg_arry  = new ( (/2,time_tot/), typeof(TVeg))
; TVeg_arry(0,:) = TVeg
; TVeg_arry(1,:) = TVeg_rc*x + TVeg_rc@yintercept
res@tiYAxisString  = f->TVeg_tavg@long_name+" (mm/0.5hr)"
plot = gsn_csm_xy(wks,x,TVeg,res)

print("plot esoil")
; ESoil_rc    = regline(x, ESoil)          ; slope
; ESoil_arry  = new ( (/2,time_tot/), typeof(ESoil))
; ESoil_arry(0,:) = ESoil
; ESoil_arry(1,:) = ESoil_rc*x + ESoil_rc@yintercept
res@tiYAxisString  = f->ESoil_tavg@long_name+" (mm/0.5hr)"
plot = gsn_csm_xy(wks,x,ESoil,res)

print("plot surface runoff")
; Qs_rc    = regline(x, Qs)          ; slope
; Qs_arry  = new ( (/2,time_tot/), typeof(Qs))
; Qs_arry(0,:) = Qs
; Qs_arry(1,:) = Qs_rc*x + Qs_rc@yintercept
res@tiYAxisString  = f->Qs_tavg@long_name+" (mm/0.5hr)"
plot = gsn_csm_xy(wks,x,Qs,res)

print("plot subsurface runoff")
; Qsb_rc    = regline(x, Qsb)          ; slope
; Qsb_arry  = new ( (/2,time_tot/), typeof(Qsb))
; Qsb_arry(0,:) = Qsb
; Qsb_arry(1,:) = Qsb_rc*x + Qsb_rc@yintercept
res@tiYAxisString  = f->Qsb_tavg@long_name+" (mm/0.5hr)"
plot = gsn_csm_xy(wks,x,Qsb,res)

print("plot aquifer recharge")
; Qrec_rc    = regline(x, Qrec)          ; slope
; Qrec_arry  = new ( (/2,time_tot/), typeof(Qrec))
; Qrec_arry(0,:) = Qrec
; Qrec_arry(1,:) = Qrec_rc*x + Qrec_rc@yintercept
res@tiYAxisString  = f->Qrec_tavg@long_name+" (mm/0.5hr)"
plot = gsn_csm_xy(wks,x,Qrec,res)


if hydr_param then
  print("plot aquifer moisture")
  ; GWwb_rc    = regline(x, GWwb)          ; slope
  ; GWwb_arry  = new ( (/2,time_tot/), typeof(GWwb))
  ; GWwb_arry(0,:) = GWwb
  ; GWwb_arry(1,:) = GWwb_rc*x + GWwb_rc@yintercept
  ; res@trYMaxF = 0.5
  ; res@trYMinF = 0.
  res@tiYAxisString  = "Aquifer Moisture (m3/m3)"
  plot = gsn_csm_xy(wks,x,GWwb,res)
  ; delete(res@trYMaxF)
  ; delete(res@trYMinF)
  ; print("plot fwsoil")
  ; ; fwsoil_rc    = regline(x, Fwsoil)          ; slope
  ; ; fwsoil_arry  = new ( (/2,time_tot/), typeof(Fwsoil))
  ; ; fwsoil_arry(0,:) = Fwsoil
  ; ; fwsoil_arry(1,:) = fwsoil_rc*x + fwsoil_rc@yintercept
  ; res@tiYAxisString  = f->fwsoil_inst@long_name+" ("+f->fwsoil_inst@units+")"
  ; plot = gsn_csm_xy(wks,x,Fwsoil,res)
  print("plot water table")
  ; WatTab_rc    = regline(x, WatTab)          ; slope
  ; WatTab_arry  = new ( (/2,time_tot/), typeof(WatTab))
  ; WatTab_arry(0,:) = WatTab
  ; WatTab_arry(1,:) = WatTab_rc*x + WatTab_rc@yintercept
  ; res@trYMaxF = 12.
  ; res@trYMinF = 0.
  res@tiYAxisString  = f->WaterTableD_tavg@long_name+" ("+f->WaterTableD_tavg@units+")"
  plot = gsn_csm_xy(wks,x,WatTab,res)
  ; delete(res@trYMaxF)
  ; delete(res@trYMinF)
else
  print("plot snow depth")
  ; SnowDepth_rc    = regline(x, SnowDepth)          ; slope
  ; SnowDepth_arry  = new ( (/2,time_tot/), typeof(SnowDepth))
  ; SnowDepth_arry(0,:) = SnowDepth
  ; SnowDepth_arry(1,:) = SnowDepth_rc*x + SnowDepth_rc@yintercept
  res@tiYAxisString  = f->SnowDepth_inst@long_name+" ("+f->SnowDepth_inst@units+")"
  plot = gsn_csm_xy(wks,x,SnowDepth,res)
end if

;
; do soil = 0,5
;   print("plot soil temp")
;   res@trYMaxF = 50.
;   res@trYMinF = 0.
;   ; SoilTemp_rc    = regline(x, SoilTemp(:,soil))          ; slope
;   ; SoilTemp_arry  = new ( (/2,time_tot/), typeof(SoilTemp))
;   ; SoilTemp_arry(0,:) = SoilTemp(:,soil)
;   ; SoilTemp_arry(1,:) = SoilTemp_rc*x + SoilTemp_rc@yintercept
;   res@tiYAxisString  = "Layer "+(soil+1)+" Soil Temp (C)"
;   plot = gsn_csm_xy(wks,x,SoilTemp(:,:,soil),res)
;   ; delete(SoilTemp_rc)
;   ; delete(SoilTemp_arry)
;   delete(res@trYMaxF)
;   delete(res@trYMinF)
; end do
do soil = 0,5
  print("plot soil moisture")
  res@trYMaxF = 0.5
  res@trYMinF = 0.
  ; SoilMoist_rc    = regline(x, SoilMoist(:,soil))          ; slope
  ; SoilMoist_arry  = new ( (/2,time_tot/), typeof(SoilMoist))
  ; SoilMoist_arry(0,:) = SoilMoist(:,soil)
  ; SoilMoist_arry(1,:) = SoilMoist_rc*x + SoilMoist_rc@yintercept
  res@tiYAxisString  = "Layer "+(soil+1)+" Soil Moisture (m3/m3)"
  plot = gsn_csm_xy(wks,x,SoilMoist(:,:,soil),res)
  ; delete(SoilMoist_rc)
  ; delete(SoilMoist_arry)
  delete(res@trYMaxF)
  delete(res@trYMinF)
end do
;
; print("plot surface temp")
; res@trYMaxF = 50.
; res@trYMinF = 0.
; ; AvgSurfT_rc    = regline(x, AvgSurfT)          ; slope
; ; AvgSurfT_arry  = new ( (/2,time_tot/), typeof(AvgSurfT))
; ; AvgSurfT_arry(0,:) = AvgSurfT
; ; AvgSurfT_arry(1,:) = AvgSurfT_rc*x + AvgSurfT_rc@yintercept
; res@tiYAxisString  = f->AvgSurfT_tavg@long_name+" (C)"
; plot = gsn_csm_xy(wks,x,AvgSurfT,res)
; delete(res@trYMaxF)
; delete(res@trYMinF)

; SoilWet_rc    = regline(x, SoilWet)          ; slope
; SoilWet_arry  = new ( (/2,time_tot/), typeof(SoilWet))
; SoilWet_arry(0,:) = SoilWet
; SoilWet_arry(1,:) = SoilWet_rc*x + SoilWet_rc@yintercept
; res@tiYAxisString = f->SoilWet_inst@long_name+" ("+f->SoilWet_inst@units+")"
; plot = gsn_csm_xy(wks,x,SoilWet_arry,res)

end
