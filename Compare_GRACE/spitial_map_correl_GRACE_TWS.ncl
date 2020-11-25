;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

; ___________________________________ Time _____________________________________
year_s           = 2000
year_e           = 2019
GRACE_start_year = 2002
GRACE_end_year   = 2019
year_sum         = year_e-year_s+1
total            = year_sum*12
Soil_thickness   = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers and 1 aquifer

; __________________________________ Path ______________________________________
path      = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"
case_name = (/"CTL","FREE_DRAIN","PUMPx10","PUMPx100"/)
case_sum  = dimsizes(case_name)
case_path = new(case_sum,"string")
row_sum   = 2
col_sum   = 2


do case_num = 0, case_sum -1
   case_path(case_num) = path+case_name(case_num)+"/outputs-spinup200yr"
end do

; ____________________________ Options _______________________________
deseason        = True
; print("deseason = " deseason)
scale           = "SE-AU" ;"AU"; "Global"; "SE-AU"
PRE             = False ; True using Precipitation-Runoff-Evaporation vs GRACE
                        ; False using SoilMoist+GWMoist+CanopInt+SWE vs GRACE
GRACE_opt       = "TWS" ;"Change_Rate"

if GRACE_opt .eq. "Change_Rate" then
  Change_Rate_opt = "CLASS" ;"Sanaa"
else if GRACE_opt .eq. "TWS" then
  output          = "daily" ; "monthly"
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


; ___________________________ Plot Check ______________________________
pic = "check"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"ViBlGrWhYeOrRe")

res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0

res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ����ҳ
res@lbLabelBarOn       = True
res@cnLinesOn          = False

;_______________________________ PFT ___________________________________
filename = case_path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f        = addfile (filename,"r")
iveg     = f->iveg({A:B},{C:D})
lat_size = dimsizes(f->y({A:B}))
lon_size = dimsizes(f->x({C:D}))
lat_crop = f->y({A:B})
lon_crop = f->x({C:D})

delete(f)
delete(filename)

;______________________________ Read Data __________________________________
if GRACE_opt .eq. "Change_Rate" then
   GRACE   = new((/total,360,720/),float)
   GRACE   = 0.0
   TWS     = new((/case_sum,total,360,720/),float)
   TWS     = 0.0
   counter = 0

   do year = year_s,year_e

     if Change_Rate_opt .eq. "CLASS" then
         if ( mod(year,4) .eq. 0) then
            dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
         else
            dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
         end if
         filename = "/g/data/w35/mm3972/data/CLASS_v1/CLASS_v1.1_"+year+".nc"
         f = addfile (filename,"r")
         do month = 0,11
            GRACE(counter,:,:) = lonFlip(f->dw(month,:,:))*3600*24*dom(month)
            counter = counter +1
         end do ;month
     else if Change_Rate_opt .eq. "Sanaa" then
         filename = "/g/data/w35/mm3972/data/GRACE/GRACE_ds_"+year+".nc"
         f = addfile (filename,"r")
         do month = 0,11
            GRACE(counter,:,:) = lonFlip(f->dw(month,:,:))*10.
            counter = counter +1
         end do ;month
      end if
      end if
   end do ;year

   GRACE@_FillValue = f->dw@_FillValue
   GRACE!0 = "time"
   GRACE!1 = "lat"
   GRACE!2 = "lon"
   GRACE&lat = f->dw&y
   GRACE&lon = f->dw&x
   delete(f)

   ;______________________________ CABLE INPUT _________________________________
   TWS_200212 = new((/case_sum,360,720/),float)
   TWS_200212 = 0.0

   do case_num = 0,case_sum-1

     print(case_name(case_num))
     counter = 0

     ; read in 2002-12
     filename = case_name(case_num)+"/cable_out_2002.nc"
     f = addfile (filename,"r")

     if PRE .eq. False then
       TWS_200212(case_num,:,:) = f->GWMoist(11,:,:)*f->GWdz*1000. + f->SWE(11,:,:) + f->CanopInt(11,:,:)
       do layer = 0,5
         TWS_200212(case_num,:,:) = TWS_200212(case_num,:,:) + f->SoilMoist(11,layer,:,:)*Soil_thickness(layer)*1000.
       end do
     end if

     ; read in 2003-2009
     do year = year_s,year_e
       print("year = " +year)
       if ( mod(year,4) .eq. 0) then
          dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
       else
          dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
       end if

       filename = case_name(case_num)+"/cable_out_"+year+".nc"
       f = addfile (filename,"r")

       do month = 0,11
         if PRE .eq. True then
           TWS(case_num,counter,:,:) = (f->Rainf(month,:,:) - f->Evap(month,:,:) - f->Qs(month,:,:) - f->Qsb(month,:,:)) *3600*24*dom(month)
         else
           TWS(case_num,counter,:,:) =  f->SWE(month,:,:) + f->CanopInt(month,:,:) + f->GWMoist(month,:,:)*f->GWdz*1000.
           do layer = 0,5
             TWS(case_num,counter,:,:) = TWS(case_num,counter,:,:) + f->SoilMoist(month,layer,:,:)*Soil_thickness(layer)*1000.
           end do
         end if
         counter = counter + 1
       end do ;month
     end do ;year

   end do ;case

   ; Calculate TWS_change
   TWS_change = new((/case_sum,total,360,720/),float)
   TWS_change = 0.0

   if PRE .eq. True then
     TWS_change = TWS
   else
     TWS_change(:,0,:,:) = TWS(:,0,:,:) - TWS_200212(:,:,:)
     TWS_change(:,1:total-1,:,:) = TWS(:,1:total-1,:,:) - TWS(:,0:total-2,:,:)
   end if
   delete(TWS)
   delete(TWS_200212)

   TWS_change!0 = "case"
   TWS_change!1 = "time"
   TWS_change!2 = "lat"
   TWS_change!3 = "lon"
   TWS_change&lat = f->Qle&y
   TWS_change&lon = f->Qle&x
   TWS_change@_FillValue   = f->Qle@_FillValue

   delete(f)
   delete(counter)

else if GRACE_opt .eq. "TWS" then

    ; ______________________________ GRACE __________________________________
    filename   = "/g/data/w35/Shared_data/Observations/GRACE/GRACE_JPLRL06M_MASCON/GRCTellus.JPL.2002_2020_gap-fill_MMY.nc"
    f          = addfile (filename, "r")

    GRACE_lat_size = dimsizes(f->lwe_thickness(0,{A:B},0))
    GRACE_lon_size = dimsizes(f->lwe_thickness(0,0,{C:D}))
    GRACE          = new((/total, GRACE_lat_size, GRACE_lon_size/), float)
    GRACE          = 0.0

    ; find time_start and time_end in GRACE input
    tot_grace        = (GRACE_end_year-GRACE_start_year+1)*12
    compensate_month = (GRACE_start_year - year_s)*12

    ; give value
    GRACE(compensate_month:compensate_month+tot_grace-1,:,:) = \
                                 f->lwe_thickness(0:tot_grace-1,{A:B},{C:D})*10.

    replace_ieeenan (GRACE, GRACE@_FillValue, 0)

    GRACE!0   = "time"
    GRACE!1   = "lat"
    GRACE!2   = "lon"
    GRACE&lat = f->lat({A:B})
    GRACE&lon = f->lon({C:D})
    delete(f)

    ;______________________________ CABLE  __________________________________
    TWS           = new((/case_sum,total,lat_size,lon_size/),float)
    TWS_ave       = new((/case_sum,lat_size,lon_size/),float)
    Rain          = new((/case_sum,total,lat_size,lon_size/),float)

    TWS           = 0.0
    TWS_ave       = 0.0
    Rain          = 0.0

    do case_num = 0,case_sum-1

       print(case_name(case_num))
       counter = 0

       do year = year_s,year_e
          print("year = " +year)

          filename = case_path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
          f        = addfile (filename,"r")

          if output .eq. "monthly" then
            do month = 0,11
               TWS(case_num,counter,:,:)       = f->GWMoist(month,{A:B},{C:D})*f->GWdz*1000. + f->CanopInt(month,{A:B},{C:D}) + f->SWE(month,{A:B},{C:D})
               Rain(case_num,counter,:,:)      = f->Rainf(month,{A:B},{C:D})
               do layer = 0,5
                  TWS(case_num,counter,:,:)    = TWS(case_num,counter,:,:) + f->SoilMoist(month,layer,{A:B},{C:D})*Soil_thickness(layer)*1000.
               end do
               counter = counter + 1
            end do ;month
          else if output .eq. "daily" then
            if ( mod(year,4) .eq. 0) then
               dom = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
            else
               dom = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
            end if

            do month = 0,11
               TWS(case_num,counter,:,:) = dim_avg_n_Wrap(f->GWMoist(dom(month):dom(month+1)-1,{A:B},{C:D}),0) \
                                          *f->GWdz*1000. + dim_avg_n_Wrap(f->CanopInt(dom(month):dom(month+1)-1,{A:B},{C:D}),0)\
                                          + dim_avg_n_Wrap(f->SWE(dom(month):dom(month+1)-1,{A:B},{C:D}), 0)
               Rain(case_num,counter,:,:)= dim_avg_n_Wrap(f->Rainf(dom(month):dom(month+1)-1,{A:B},{C:D}), 0)
               do layer = 0,5
                  TWS(case_num,counter,:,:) = TWS(case_num,counter,:,:) + dim_avg_n_Wrap( \
                                             f->SoilMoist(dom(month):dom(month+1)-1,layer,{A:B},{C:D}),0)\
                                             *Soil_thickness(layer)*1000.
               end do
               counter = counter + 1
            end do ;month
          end if ; output
          end if ; output
       end do ;year
    end do ;case

    ; Calculate 2004-2009 TWS average
    tot_TWS       = (2009-2004+1) * 12 -1
    TWS_ts_str_yr = (2004 - year_s) * 1
    TWS_ave =  dim_avg_n_Wrap(TWS(:, TWS_ts_str_yr:TWS_ts_str_yr + tot_TWS,:,:),1)

    ; Calculate TWS anormaly
    do time = 0,total -1
      TWS(:,time,:,:)      = TWS(:,time,:,:) - TWS_ave
    end do

    TWS_regid = linint2(GRACE&lon,GRACE&lat,TWS,True,lon_crop,lat_crop,0)

    res@cnLevels = fspan(-100,100,10)
    plot  = gsn_csm_contour(wks,dim_avg_n_Wrap(GRACE(:,:,:),0),res)
    plot  = gsn_csm_contour(wks,dim_avg_n_Wrap(TWS_regid(0,:,:,:),0),res)
    plot  = gsn_csm_contour(wks,dim_avg_n_Wrap(TWS(0,:,:,:),0),res)
    GRACE = where( .not. ismissing(TWS_regid(0,:,:,:)), GRACE, GRACE@_FillValue)

    delete(TWS)
    delete(TWS_ave)
    delete(res@cnLevels)

    TWS_regid!0   = "case"
    TWS_regid!1   = "time"
    TWS_regid!2   = "lat"
    TWS_regid!3   = "lon"
    TWS_regid&lat = GRACE&lat
    TWS_regid&lon = GRACE&lon
    TWS_regid@_FillValue = f->Qle@_FillValue

    delete(counter)
    delete(f)

end if ; GRACE_opt
end if ; GRACE_opt


; __________________________ Deseasonalize _____________________________
; rmMonAnnCycTLL : Removes the annual cycle from "monthly" data.
if deseason .eq. True then
  print(deseason)

  GRACE = rmMonAnnCycTLL(GRACE) ; dimensioned time x lat x lon

  do case_num = 0,case_sum-1
    if GRACE_opt .eq. "Change_Rate" then
       TWS_change(case_num,:,:,:) = rmMonAnnCycTLL(TWS_change(case_num,:,:,:))
    else if GRACE_opt .eq. "TWS" then
       TWS_regrid(case_num,:,:,:) = rmMonAnnCycTLL(TWS_regrid(case_num,:,:,:))
    end if
    end if
  end do
end if

plot  = gsn_csm_contour(wks,dim_avg_n_Wrap(GRACE(:,:,:),0),res)
plot  = gsn_csm_contour(wks,dim_avg_n_Wrap(TWS_regid(0,:,:,:),0),res)

; _____________________  CALCULATE CORRELATION ________________________
GRACE_TWS  = new((/case_sum,GRACE_lat_size,GRACE_lon_size/),float)
GRACE_TWS  = 0.0

do case_num = 0,case_sum-1

  if GRACE_opt .eq. "Change_Rate" then
      GRACE_TWS(case_num,:,:) = escorc(TWS_change(case|case_num,lat|:,lon|:,time|:), GRACE(lat|:,lon|:,time|:))
      delete(TWS_change)
  else if GRACE_opt .eq. "TWS" then
      GRACE_TWS(case_num,:,:) = escorc(TWS_regrid(case|case_num,lat|:,lon|:,time|:), GRACE(lat|:,lon|:,time|:))
      delete(TWS)
  end if
  end if
; escorc: Computes the (Pearson) sample linear cross-correlations at lag 0 only.
; The rightmost dimension is usually time.
end do

GRACE_TWS!0 = "case"
GRACE_TWS!1 = "lat"
GRACE_TWS!2 = "lon"
GRACE_TWS&lat = GRACE&lat
GRACE_TWS&lon = GRACE&lon

plot = gsn_csm_contour(wks,GRACE_TWS(0,:,:)*1.0,res)
delete(GRACE)
delete(TWS_regrid)
delete(res)
delete(pic)
delete(wks)

;_____________________________ PLOT ___________________________________
if deseason .eq. True then
   pic = "./plots/"+GRACE_opt+"-GRACE_TWS_correl_"+scale+"_deseasonal_"+year_s+"-"+year_e
else
   pic = "./plots/"+GRACE_opt+"-GRACE_TWS_correl_"+scale+"_"+year_s+"-"+year_e
end if

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

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
  res@tmYLLabels      = (/"40~S~o~N~N","35~S~o~N~N","30~S~o~N~","25~S~o~N~S"/)
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = (/140,145,150/)
  res@tmXBLabels      = (/"140~S~o~N~W","145~S~o~N~W","150~S~o~N~"/)
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
pres                    = True                                      ; ���� panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.013
pres@lbLabelFontHeightF = 0.013
pres@lbAutoManage       = False
pres@lbTitleDirection   = "Across"             ; title direction
pres@lbTitlePosition    = "Bottom";Right"
;pres@pmLabelBarSide     = "Bottom";"Right"
;pres@pmLabelBarOrthogonalPosF = -0.021
;pres@pmLabelBarParallelPosF = 0.

plots = new(case_sum,graphic)

; ____________________ all PFT __________________________
print("for all PFTs")
do case_num = 0,case_sum-1
    print(case_name(case_num)+" "+avg(GRACE_TWS(case_num,{A:B},{C:D})))
    res@tiMainString = case_name(case_num)
    res@cnLevels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
    plots(case_num) = gsn_csm_contour_map(wks,GRACE_TWS(case_num,{A:B},{C:D}),res)
end do
gsn_panel(wks,(/plots/),(/row_sum,col_sum/),pres)

;
; ; ____________________ forest __________________________
; ; tmp var
; tmp   = new((/case_sum, GRACE_lat_size, GRACE_lon_size/),float)
; tmp   = 0.0
; tmp!0 = "case"
; tmp!1 = "lat"
; tmp!2 = "lon"
; tmp&lat = GRACE_TWS&lat
; tmp&lon = GRACE_TWS&lon
; tmp@_FillValue = GRACE_TWS@_FillValue
; print("for forest")
; pres@txString = "forest"
; do case_num = 0,case_sum-1
;   tmp(case_num,:,:) = where(iveg .ge. 1 .and. iveg .le. 4, GRACE_TWS(case_num,:,:), GRACE_TWS@_FillValue )
;   res@tiMainString = case_name(case_num)
;   if scale .eq. "Global" then
;     print(case_name_brf(case_num)+" "+avg(tmp(case_num,:,:)))
;     levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
;     res@cnLevels = levels
;     plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,:,:),res)
;   else if scale .eq. "AU" then
;     print(case_name_brf(case_num)+" "+avg(tmp(case_num,{-40:-27},{140:155})))
;     levels = (/0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6/)
;     res@cnLevels = levels
;     plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,{-40:-27},{140:155}),res)
;   end if
;   end if
; end do
; gsn_panel(wks,(/plots/),(/1,1/),pres)
;
;
; ; ____________________ grass+shrub+crop __________________________
; print("for grass+shrub+crop")
; pres@txString = "grass+shrub+crop"
; do case_num = 0,case_sum-1
;   tmp(case_num,:,:) = where(iveg .ge. 5 .and. iveg .le. 11, GRACE_TWS(case_num,:,:), GRACE_TWS@_FillValue )
;   res@tiMainString = case_name_brf(case_num)
;   if scale .eq. "Global" then
;     print(case_name_brf(case_num)+" "+avg(tmp(case_num,:,:)))
;     levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
;     res@cnLevels = levels
;     plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,:,:),res)
;   else if scale .eq. "AU" then
;     print(case_name_brf(case_num)+" "+avg(tmp(case_num,{-40:-27},{140:155})))
;     levels = (/0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6/)
;     res@cnLevels = levels
;     plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,{-40:-27},{140:155}),res)
;   end if
;   end if
; end do
; gsn_panel(wks,(/plots/),(/1,1/),pres)

end
