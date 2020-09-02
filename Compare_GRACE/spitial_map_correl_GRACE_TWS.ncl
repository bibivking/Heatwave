;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3/outputs/GSWP/"\
             /)

case_name_brf = (/\
              "GSWP"\
               /)

var_name= (/\
           "SoilMoist",\     ;"m^3/m^3"  Average layer soil moisture
           "GWMoist",\       ;"mm3/mm3" Aquifer mositure content, Groundwater_off doesn't have
           "CanopInt",\      ;"kg/m^2", "Canopy intercepted water storage"
           "SWE",\           ;"kg/m^2", "Snow water equivalent"
           "GWdz"\           ;"Mean aquifer layer thickness" ;
         	 /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers and 1 aquifer

; ____________________________ Options _______________________________
deseason = True
print(deseason)
scale = "World" ;"AU" or "World"
PRE   = False ; True using Precipitation-Runoff-Evaporation vs GRACE
             ; False using SoilMoist+GWMoist+CanopInt+SWE vs GRACE
GRACE_opt       = "TWS";"Change_Rate"

if GRACE_opt .eq. "Change_Rate" then
  Change_Rate_opt = "CLASS" ;"Sanaa"
else if GRACE_opt .eq. "TWS" then
  output          = "daily"
  use_month       = "May"
end if
end if


;______________________________ Data _________________________________
year_s = 2003
year_e = 2010
year_sum = year_e-year_s+1
total    = year_sum*12
var_sum  = dimsizes(var_name)
case_sum = dimsizes(case_name)

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
filename = case_name(0)+"/cable_out_"+year_s+".nc"
f = addfile (filename,"r")
iveg = f->iveg
delete(f)
delete(filename)

;______________________________ Read Data __________________________________
GRACE      = new((/total,360,720/),float)
GRACE      = 0.0

if GRACE_opt .eq. "Change_Rate" then

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

    filename   = "/g/data/w35/mm3972/data/GRACE/GRACE_JPL_RL06/GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc"
    f          = addfile (filename,"r")

    GRACE(0:4,:,:)       = dble2flt(f->lwe_thickness(7:11,:,:))*10.
    GRACE(6:total-1,:,:) = dble2flt(f->lwe_thickness(12:total+5,:,:))*10.
    GRACE@_FillValue     = dble2flt(f->lwe_thickness@_FillValue)

    May_ave = new((/360,720/),float)
    Jun_ave = new((/360,720/),float)
    Jul_ave = new((/360,720/),float)
    May_ave = 0.
    Jun_ave = 0.
    Jul_ave = 0.

    do year = year_s+1,year_e
      print("add year :" + year)
      May_ave = May_ave + GRACE((year-year_s)*12+4,:,:)
      Jun_ave = Jun_ave + GRACE((year-year_s)*12+5,:,:)
      Jul_ave = Jul_ave + GRACE((year-year_s)*12+6,:,:)
    end do

    May_ave = May_ave/(year_e-year_s)
    Jun_ave = Jun_ave/(year_e-year_s)
    Jul_ave = Jul_ave/(year_e-year_s)

    if use_month .eq. "Jul" then
       GRACE(5,:,:) = GRACE(6,:,:) * (Jun_ave/Jul_ave)
    else if use_month .eq. "May" then
       GRACE(5,:,:) = GRACE(4,:,:) * (Jun_ave/May_ave)
    end if
    end if

    delete(May_ave)
    delete(Jun_ave)
    delete(Jul_ave)

    GRACE!0 = "time"
    GRACE!1 = "lat"
    GRACE!2 = "lon"
    GRACE&lat = dble2flt(f->lat)
    GRACE&lon = dble2flt(f->lon)
    delete(f)

    ;______________________________ CABLE INPUT __________________________________
    TWS           = new((/case_sum,total,360,720/),float)
    TWS_ave       = new((/case_sum,360,720/),float)
    Rain          = new((/case_sum,total,360,720/),float)

    TWS           = 0.0
    TWS_ave       = 0.0
    Rain          = 0.0

    do case_num = 0,case_sum-1

       print(case_name(case_num))
       counter = 0

       do year = year_s,year_e
          print("year = " +year)

          filename = case_name(case_num)+"/cable_out_"+year+".nc"
          f = addfile (filename,"r")

          if output .eq. "monthly" then

            do month = 0,11
               TWS(case_num,counter,:,:)       = f->GWMoist(month,:,:)*f->GWdz*1000. + f->CanopInt(month,:,:) + f->SWE(month,:,:)
               Rain(case_num,counter,:,:)      = f->Rainf(month,:,:)
               do layer = 0,5
                  TWS(case_num,counter,:,:)    = TWS(case_num,counter,:,:) + f->SoilMoist(month,layer,:,:)*Soil_thickness(layer)*1000.
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
               TWS(case_num,counter,:,:) = dim_avg_n_Wrap(f->GWMoist(dom(month):dom(month+1)-1,:,:),0) \
                                          *f->GWdz*1000. + dim_avg_n_Wrap(f->CanopInt(dom(month):dom(month+1)-1,:,:),0)\
                                          + dim_avg_n_Wrap(f->SWE(dom(month):dom(month+1)-1,:,:), 0)
               Rain(case_num,counter,:,:)= dim_avg_n_Wrap(f->Rainf(dom(month):dom(month+1)-1,:,:), 0)
               do layer = 0,5
                  TWS(case_num,counter,:,:) = TWS(case_num,counter,:,:) + dim_avg_n_Wrap( \
                                             f->SoilMoist(dom(month):dom(month+1)-1,layer,:,:),0)\
                                             *Soil_thickness(layer)*1000.
               end do
               counter = counter + 1
            end do ;month
          end if
          end if
       end do ;year
    end do ;case

    ; Calculate 2004-2009 TWS average
    GRACE_end_year = 2009
    if year_e .gt. 2009 then
      adjust  = (year_e - GRACE_end_year)*12
    else
      adjust  = 0
      GRACE_end_year = year_e
    end if
    TWS_ave        =  dim_avg_n_Wrap(TWS(:,12:total-1 - adjust,:,:),1)

    ; Calculate TWS anormaly
    do time = 0,total -1
      TWS(:,time,:,:)      = TWS(:,time,:,:) - TWS_ave
    end do

    delete(TWS_ave)
    levels = fspan(-100,100,10)
    res@cnLevels = levels
    plot  = gsn_csm_contour(wks,GRACE(0,:,:)*1.0,res)
    GRACE = where( .not. ismissing(TWS(0,:,:,:)), GRACE, GRACE@_FillValue)
    delete(levels)
    delete(res@cnLevels)

    TWS!0 = "case"
    TWS!1 = "time"
    TWS!2 = "lat"
    TWS!3 = "lon"
    TWS&lat = f->Qle&y
    TWS&lon = f->Qle&x
    TWS@_FillValue   = f->Qle@_FillValue

    delete(counter)
    delete(f)

end if
end if

plot = gsn_csm_contour(wks,TWS(0,0,:,:)*1.0,res)

; __________________________ Deseasonalize _____________________________
; rmMonAnnCycTLL : Removes the annual cycle from "monthly" data.
if deseason .eq. True then
  print(deseason)

  GRACE = rmMonAnnCycTLL(GRACE) ; dimensioned time x lat x lon

  do case_num = 0,case_sum-1
    if GRACE_opt .eq. "Change_Rate" then
       TWS_change(case_num,:,:,:) = rmMonAnnCycTLL(TWS_change(case_num,:,:,:))
    else if GRACE_opt .eq. "TWS" then
       TWS(case_num,:,:,:) = rmMonAnnCycTLL(TWS(case_num,:,:,:))
    end if
    end if
  end do
end if

plot = gsn_csm_contour(wks,GRACE(0,:,:)*1.0,res)
plot = gsn_csm_contour(wks,TWS(0,0,:,:)*1.0,res)

; _____________________  CALCULATE CORRELATION ________________________
GRACE_TWS  = new((/case_sum,360,720/),float)
GRACE_TWS  = 0.0

do case_num = 0,case_sum-1

  if GRACE_opt .eq. "Change_Rate" then
      GRACE_TWS(case_num,:,:) = escorc(TWS_change(case|case_num,lat|:,lon|:,time|:), GRACE(lat|:,lon|:,time|:))
      delete(TWS_change)
  else if GRACE_opt .eq. "TWS" then
      GRACE_TWS(case_num,:,:) = escorc(TWS(case|case_num,lat|:,lon|:,time|:), GRACE(lat|:,lon|:,time|:))
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

if scale .eq. "World" then
   res@trYMinF            = -60
   res@trYMaxF            = 80
   res@mpMinLatF          = -60
   res@mpMaxLatF          = 80
   res@mpCenterLonF       = 160

   res@tmYLMode  = "Explicit"             ;��Y����������������,��ʾ�Լ���Ҫ��ֵ
   res@tmYLValues=(/60,30,0,-30,-60/)
   res@tmYLLabels=(/"60~S~o~N~N","30~S~o~N~N","0~S~o~N~","30~S~o~N~S","60~S~o~N~S"/)
   res@tmXBMode  = "Explicit"
   res@tmXBValues=(/240,300,0,60,120,180/)
   res@tmXBLabels=(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)

else if scale .eq. "AU" then

   res@gsnAddCyclic       = False
   res@trYMinF            = -40
   res@trYMaxF            = -27
   res@mpMinLatF          = -40
   res@mpMaxLatF          = -27
   res@trXMinF            = 140
   res@trXMaxF            = 155
   res@mpMinLonF          = 140
   res@mpMaxLonF          = 155

   res@tmYLMode  = "Explicit"             ;��Y����������������,��ʾ�Լ���Ҫ��ֵ
   res@tmYLValues=(/-30,-40/)
   res@tmYLLabels=(/"30~S~o~N~S","40~S~o~N~S"/)
   res@tmXBMode  = "Explicit"
   res@tmXBValues=(/140,150/)
   res@tmXBLabels=(/"140~S~o~N~E","150~S~o~N~E"/)
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

res@tiMainString = scale+ " "+case_name_brf(case_num)

if scale .eq. "World" then
   print(case_name_brf(case_num)+" "+avg(GRACE_TWS(case_num,:,:)))
   levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
   res@cnLevels = levels
   plots(case_num) = gsn_csm_contour_map(wks,GRACE_TWS(case_num,:,:),res)
else if scale .eq. "AU" then
   print(case_name_brf(case_num)+" "+avg(GRACE_TWS(case_num,{-40:-27},{140:155})))
   levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
   res@cnLevels = levels
   plots(case_num) = gsn_csm_contour_map(wks,GRACE_TWS(case_num,{-40:-27},{140:155}),res)
end if
end if

end do
gsn_panel(wks,(/plots/),(/1,1/),pres)

; tmp var
tmp   = new((/case_sum,360,720/),float)
tmp   = 0.0
tmp!0 = "case"
tmp!1 = "lat"
tmp!2 = "lon"
tmp&lat = GRACE_TWS&lat
tmp&lon = GRACE_TWS&lon
tmp@_FillValue = GRACE_TWS@_FillValue

; ____________________ forest __________________________
print("for forest")
pres@txString = "forest"
do case_num = 0,case_sum-1
  tmp(case_num,:,:) = where(iveg .ge. 1 .and. iveg .le. 4, GRACE_TWS(case_num,:,:), GRACE_TWS@_FillValue )
  res@tiMainString = case_name_brf(case_num)
  if scale .eq. "World" then
    print(case_name_brf(case_num)+" "+avg(tmp(case_num,:,:)))
    levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
    res@cnLevels = levels
    plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,:,:),res)
  else if scale .eq. "AU" then
    print(case_name_brf(case_num)+" "+avg(tmp(case_num,{-40:-27},{140:155})))
    levels = (/0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6/)
    res@cnLevels = levels
    plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,{-40:-27},{140:155}),res)
  end if
  end if
end do
gsn_panel(wks,(/plots/),(/1,1/),pres)


; ____________________ grass+shrub+crop __________________________
print("for grass+shrub+crop")
pres@txString = "grass+shrub+crop"
do case_num = 0,case_sum-1
  tmp(case_num,:,:) = where(iveg .ge. 5 .and. iveg .le. 11, GRACE_TWS(case_num,:,:), GRACE_TWS@_FillValue )
  res@tiMainString = case_name_brf(case_num)
  if scale .eq. "World" then
    print(case_name_brf(case_num)+" "+avg(tmp(case_num,:,:)))
    levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
    res@cnLevels = levels
    plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,:,:),res)
  else if scale .eq. "AU" then
    print(case_name_brf(case_num)+" "+avg(tmp(case_num,{-40:-27},{140:155})))
    levels = (/0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6/)
    res@cnLevels = levels
    plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,{-40:-27},{140:155}),res)
  end if
  end if
end do
gsn_panel(wks,(/plots/),(/1,1/),pres)

end
