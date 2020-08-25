;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 1000000000 ;
end setvalues
;_____________________________________

begin

case_name = (/\
              "/g/data/w35/mm3972/model/cable/runs/my_version/run_GSWP3/outputs/GSWP/"\
             /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

; ______________________________ Pre-define ________________________
year_s   = 2003
year_e   = 2008
year_sum = year_e-year_s+1
total    = year_sum*12
case_sum = dimsizes(case_name)

; ____________________________ Options _______________________________

scale        = "AU" ;"Global" ; "50S-N" ; "SE-AU" ; "AU"
use_month    = "May"

if scale .eq. "Global" then
  A = -90.
  B = 60.
  C = 0.
  D = 360.
else if scale .eq. "50S-N" then
  A = -50.
  B = 50.
  C = 0.
  D = 360.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -27.
  C = 140.
  D = 155.
else if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
end if
end if
end if
end if

; _______________________ read GRACE and fill gap ____________________________
GRACE      = new((/total,360,720/),float)
uncertainty= new((/total,360,720/),float)
GRACE      = 0.
uncertainty= 0.

filename   = "/g/data/w35/mm3972/data/GRACE/GRACE_JPL_RL06/GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc"
f          = addfile (filename,"r")

GRACE(0:4,:,:)       = dble2flt(f->lwe_thickness(7:11,:,:))*10.
GRACE(6:total-1,:,:)   = dble2flt(f->lwe_thickness(12:total+5,:,:))*10.

uncertainty(0:4,:,:) = dble2flt(f->uncertainty(7:11,:,:))*10.
uncertainty(6:total-1,:,:)= dble2flt(f->uncertainty(12:total+5,:,:))*10.

GRACE@_FillValue       = dble2flt(f->lwe_thickness@_FillValue)
uncertainty@_FillValue = dble2flt(f->lwe_thickness@_FillValue)

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

;______________________________ INPUT __________________________________
TWS     = new((/case_sum,total,360,720/),float)
TWS_ave = new((/case_sum,360,720/),float)
TWS     = 0.0
TWS_ave = 0.0

do case_num = 0,case_sum-1
   print(case_name(case_num))
   counter = 0
   do year = year_s,year_e
      print("year = " +year)

      filename = case_name(case_num)+"/cable_out_"+year+".nc"
      f = addfile (filename,"r")

      do month = 0,11
         TWS(case_num,counter,:,:) = f->GWMoist(month,:,:)*f->GWdz*1000. + f->CanopInt(month,:,:) + f->SWE(month,:,:)
         do layer = 0,5
            TWS(case_num,counter,:,:) = TWS(case_num,counter,:,:) + f->SoilMoist(month,layer,:,:)*Soil_thickness(layer)*1000.
         end do
         counter = counter + 1
      end do ;month
   end do ;year
end do ;case

; Calculate 2004-2008 TWS average
TWS_ave =  dim_avg_n_Wrap(TWS(:,12:71,:,:),1)

; Calculate TWS anormaly
do time = 0,total-1
  TWS(:,time,:,:) = TWS(:,time,:,:) - TWS_ave
end do

delete(TWS_ave)

print("point 1 ")

; ____________________________ Modify GRACE data _______________________________
GRACE       = where( .not. ismissing(TWS(0,:,:,:)), GRACE, GRACE@_FillValue)
uncertainty = where( .not. ismissing(TWS(0,:,:,:)), uncertainty, uncertainty@_FillValue)

lat = f->Qle&y
lon = f->Qle&x

TWS!0 = "case"
TWS!1 = "time"
TWS!2 = "lat"
TWS!3 = "lon"
TWS&lat = lat
TWS&lon = lon
TWS@_FillValue   = f->Qle@_FillValue

GRACE!0 = "time"
GRACE!1 = "lat"
GRACE!2 = "lon"
GRACE&lat = lat
GRACE&lon = lon

uncertainty!0 = "time"
uncertainty!1 = "lat"
uncertainty!2 = "lon"
uncertainty&lat = lat
uncertainty&lon = lon

; _______________________________ Time_serial __________________________________
GRACE_time = new(total,float)
TWS_time   = new((/case_sum,total/),float)
GRACE_time = 0.
TWS_time   = 0

do time = 0, total-1
  do case_num = 0,case_sum-1
    TWS_time(case_num,time) = avg(TWS(case_num,time,{A:B},{C:D}))
  end do
  GRACE_time(time) = avg(GRACE(time,{A:B},{C:D}))
end do

;_____________________________ PLOT ___________________________________
pic = "Time_serial_GRACE_vs_CABLE_GRACE06_"+scale

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"BlueYellowRed")

res                    = True
res@gsnMaximize        = True
res@vpHeightF          = 0.6
res@vpWidthF           = 0.6
res@vpXF               = 0.3
res@vpYF               = 0.85
res@pmLegendWidthF     = 0.14                     ; Resize legend width
res@pmLegendHeightF    = 0.11                     ; and height

res@xyLineThicknesses  = (/  1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0/)          ; make second line thicker
;res@pmLegendDisplayMode = "Always"
res@xyLineColors       = (/"blue","red","green","orange","red","red","red"/)          ; change line color
res@xyExplicitLegendLabels = (/"GRACE","CABLE"/)
;res@trYMinF = 0
res@pmLegendDisplayMode = "Always"

x = ispan(1,total,1)

do case_num = 0,case_sum-1
   res@tiMainString = "GRACE vs CABLE-GSWP"
   plot = gsn_csm_xy(wks,x,(/GRACE_time, TWS_time(case_num,:)/),res)
end do

end
