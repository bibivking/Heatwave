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
case_name = "Princeton_ctl_watmove_new_veg"
path1     = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output"

; case_name = ["ERAI_ctl","ERAI_beta_exp","ERAI_gw_off","ERAI_or","ERAI_hyre","ERAI_watmove"]

; ================================= Predefine ==================================
year_s = 2008
year_e = 2008
year_sum = year_e-year_s+1
total = year_sum*366

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Soil_thickness1 = Soil_thickness
;(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE


; ================================== OPTIONS ===================================
plot_type = "all"
force_name= "Princeton";"ERAI"
time_scale= "daily" ; "annual"

; ================================ Set Regrid ==================================
interp_method = "neareststod"  ;"neareststod"; "conserve" "patch";"bilinear"

filename  = path1+"/LIS.CABLE."+year_s+"090100.d01.nc"
f = addfile (filename,"r")

lat2d    = f->lat(0,:,:)
lon2d    = f->lon(0,:,:)
lon2d    = where(lon2d.lt.0,360+lon2d,lon2d)

minlon  = ceil(min(lon2d))+0.5 ; 116.5516
minlat  = ceil(min(lat2d))+0.5 ; -44.43771
maxlon  = floor(max(lon2d))-0.5 ; 154.5981
maxlat  = floor(max(lat2d))-0.5 ; -12.5853

; study region
A = minlat
B = maxlat
C = minlon
D = maxlon

; ================================ READ CABLE ==================================
filename = path1+"/LIS.CABLE."+year_s+"010100.d01.nc"
f        = addfile (filename,"r")

if time_scale .eq. "annual" then
    EVAPFBL     = dim_sum_n_Wrap(f->EVAPFBL_tavg,0)
    WATMOVE     = dim_sum_n_Wrap(f->watmove_tavg,0)
    TVeg        = dim_sum_n_Wrap(f->TVeg_tavg,0)
    WaterTable  = dim_sum_n_Wrap(f->WaterTableD_tavg,0)/1000.
else if time_scale .eq. "daily" then
    SoilMoist   = new((/total,6, dimsizes(f->TVeg_tavg(0,:,0)),dimsizes(f->TVeg_tavg(0,0,:))/),float)
    EVAPFBL     = new((/total,6, dimsizes(f->TVeg_tavg(0,:,0)),dimsizes(f->TVeg_tavg(0,0,:))/),float)
    WATMOVE     = new((/total,6, dimsizes(f->TVeg_tavg(0,:,0)),dimsizes(f->TVeg_tavg(0,0,:))/),float)
    TVeg        = new((/total,dimsizes(f->TVeg_tavg(0,:,0)),dimsizes(f->TVeg_tavg(0,0,:))/),float)
    WaterTable  = new((/total,dimsizes(f->TVeg_tavg(0,:,0)),dimsizes(f->TVeg_tavg(0,0,:))/),float)
    GWwb        = new((/total,dimsizes(f->TVeg_tavg(0,:,0)),dimsizes(f->TVeg_tavg(0,0,:))/),float)
    Qrec        = new((/total,dimsizes(f->TVeg_tavg(0,:,0)),dimsizes(f->TVeg_tavg(0,0,:))/),float)
end if
end if

SoilMoist   = 0.
EVAPFBL     = 0.
WATMOVE     = 0.
TVeg        = 0.

delete(f)

counter = 0
do year = year_s,year_e
    print("year = " +year)

    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    ; _________________________ read LIS-CABLE ____________________________
    do month = 1, 12
      if month .lt. 10 then
        filename  = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename  = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if
      f  = addfile (filename,"r")
      print("month = " +month)
      if time_scale .eq. "annual" then
          EVAPFBL     = EVAPFBL   + dim_sum_n_Wrap(f->EVAPFBL_tavg,0)*3600*24
          WATMOVE     = WATMOVE   + dim_sum_n_Wrap(f->watmove_tavg,0)*3600*24
          TVeg        = TVeg      + dim_sum_n_Wrap(f->TVeg_tavg,0)*3600*24
      else if time_scale .eq. "daily" then
          SoilMoist(counter:counter+dom(month-1)-1,:,:,:)   = f->SoilMoist_tavg
          EVAPFBL(counter:counter+dom(month-1)-1,:,:,:)     = f->EVAPFBL_tavg*3600*24
          WATMOVE(counter:counter+dom(month-1)-1,:,:,:)     = f->watmove_tavg*3600*24
          TVeg(counter:counter+dom(month-1)-1,:,:)          = f->TVeg_tavg*3600*24
          WaterTable(counter:counter+dom(month-1)-1,:,:)    = f->WaterTableD_tavg/1000.
          GWwb(counter:counter+dom(month-1)-1,:,:)          = f->GWwb_tavg
          Qrec(counter:counter+dom(month-1)-1,:,:)          = f->Qrec_tavg*3600*24
          counter = counter+dom(month-1)
      end if
      end if
      delete(f)
    end do ; month
end do ; year

if time_scale .eq. "annual" then
    EVAPFBL     = EVAPFBL/year_sum
    WATMOVE     = WATMOVE/year_sum
    TVeg        = TVeg/year_sum
end if

;print(TVeg)
; ================================= Plotting ===================================
res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0

res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ
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

;**************** ͼ�������� *****************
pres                    = True                                      ; ���� panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
;pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01

res@gsnLeftString = ""
res@gsnRightString = ""

; ============================= PLOTTING =============================
pic1 = "./plots/LIS_water_movement_"+case_name+"_"+time_scale+"_"+year_s+"-"+year_e
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"WhiteBlueGreenYellowRed") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

if time_scale .eq. "annual" then
    plots = new(2,graphic)
    levels = fspan(0,800,41)
    res@cnLevels = levels
    res@tiMainString = "TVeg"
    plots(0)  = gsn_csm_contour(wks1,TVeg,res)
    delete(levels)
    delete(res@cnLevels)

    res@tiMainString = "WaterTable"
    plots(1)  = gsn_csm_contour(wks1,WaterTable,res)
    gsn_panel(wks1,(/plots/),(/1,2/),pres)
    delete(plots)

    ; EVAPFBL
    plots = new(6,graphic)
    pres@txString = "Transpiration from ith layer"
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    levels = fspan(0,400,21)
    res@cnLevels = levels
    do i = 0,5
      res@tiMainString = i+1 + " layer"
      plots(i) = gsn_csm_contour(wks1,EVAPFBL(i,:,:),res)
    end do
    gsn_panel(wks1,(/plots/),(/2,3/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; EVAPFBL/TVeg
    plots = new(6,graphic)
    pres@txString = "Percentage of Transpiration from ith layer"
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    levels = fspan(0,50,21)
    res@cnLevels = levels

    tmp = TVeg
    tmp = 0.

    do i = 0,5
      res@tiMainString = i+1 + " layer"
      TVeg = where(TVeg .gt.0, TVeg, -99999.)
      tmp = tmp + EVAPFBL(i,:,:)/TVeg*100.
      plots(i) = gsn_csm_contour(wks1,EVAPFBL(i,:,:)/TVeg*100.,res)
    end do
    gsn_panel(wks1,(/plots/),(/2,3/),pres)
    delete(levels)
    delete(res@cnLevels)

    ; WATMOVE
    pres@txString = "Water move from upper layer to ith layer"
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = fspan(-400,400,41)
    res@cnLevels = levels
    do i = 0,5
      res@tiMainString = i+1 + " layer"
      plots(i) = gsn_csm_contour(wks1,WATMOVE(i,:,:),res)
    end do
    gsn_panel(wks1,(/plots/),(/2,3/),pres)
    delete(levels)
    delete(res@cnLevels)
else if time_scale .eq. "daily" then

    plots = new(8,graphic)

    do time = 0, total -1

      ; levels = fspan(0,4,41)
      ; res@cnLevels = levels
      ; res@tiMainString = "TVeg"
      ; plots(0)  = gsn_csm_contour(wks1,TVeg(time,:,:),res)
      ; delete(levels)
      ; delete(res@cnLevels)
      ;
      ; ; EVAPFBL
      ; pres@txString = "Transpiration from ith layer of day-"+time
      ; res@cnFillPalette = "WhiteBlueGreenYellowRed"
      ; levels = fspan(0,2,21)
      ; res@cnLevels = levels
      ; do i = 0,5
      ;   res@tiMainString = i+1 + " layer"
      ;   plots(i+1) = gsn_csm_contour(wks1,EVAPFBL(time,i,:,:),res)
      ; end do
      ; gsn_panel(wks1,(/plots/),(/3,3/),pres)
      ; delete(levels)
      ; delete(res@cnLevels)

      res@tiMainString = "WaterTable"
      res@cnFillPalette = "WhiteBlueGreenYellowRed"
      levels = fspan(0,15,31)
      res@cnLevels = levels
      plots(0)  = gsn_csm_contour(wks1,WaterTable(time,:,:),res)
      delete(levels)
      delete(res@cnLevels)

      res@tiMainString = "Recharge to aquifer"
      res@cnFillPalette = "ViBlGrWhYeOrRe"
      levels = fspan(-0.5,0.5,21)
      res@cnLevels = levels
      plots(1)  = gsn_csm_contour(wks1,Qrec(time,:,:),res)
      delete(levels)
      delete(res@cnLevels)

      ; WATMOVE
      pres@txString = "Water move from upper layer to ith layer of day-"+time
      res@cnFillPalette = "ViBlGrWhYeOrRe"
      levels = fspan(-2,2,41)
      res@cnLevels = levels
      do i = 0,5
        res@tiMainString = i+1 + " layer"
        plots(i+2) = gsn_csm_contour(wks1,WATMOVE(time,i,:,:),res)
      end do
      gsn_panel(wks1,(/plots/),(/3,3/),pres)
      delete(levels)
      delete(res@cnLevels)

      ;
      ; res@tiMainString = "Groundwater Moisture"
      ; res@cnFillPalette = "WhiteBlueGreenYellowRed"
      ; levels = fspan(0,0.5,21)
      ; res@cnLevels = levels
      ; plots(1)  = gsn_csm_contour(wks1,GWwb(time,:,:),res)
      ; delete(levels)
      ; delete(res@cnLevels)
      ;
      ; ; SoilMoist
      ; pres@txString = "Soil Moisture day-"+time
      ; res@cnFillPalette = "WhiteBlueGreenYellowRed"
      ; levels = fspan(0,0.5,21)
      ; res@cnLevels = levels
      ; do i = 0,5
      ;   res@tiMainString = i+1 + " layer"
      ;   plots(i+2) = gsn_csm_contour(wks1,SoilMoist(time,i,:,:),res)
      ; end do
      ; gsn_panel(wks1,(/plots/),(/3,3/),pres)
      ; delete(levels)
      ; delete(res@cnLevels)


      ; ; EVAPFBL/TVeg
      ; plots = new(6,graphic)
      ; pres@txString = "Percentage of Transpiration from ith layer"
      ; res@cnFillPalette = "WhiteBlueGreenYellowRed"
      ; levels = fspan(0,50,21)
      ; res@cnLevels = levels
      ;
      ; tmp = TVeg
      ; tmp = 0.
      ;
      ; do i = 0,5
      ;   res@tiMainString = i+1 + " layer"
      ;   TVeg = where(TVeg .gt.0, TVeg, -99999.)
      ;   tmp = tmp + EVAPFBL(i,:,:)/TVeg*100.
      ;   plots(i) = gsn_csm_contour(wks1,EVAPFBL(i,:,:)/TVeg*100.,res)
      ; end do
      ; gsn_panel(wks1,(/plots/),(/2,3/),pres)
      ; delete(levels)
      ; delete(res@cnLevels)
    end do
end if
end if
end
