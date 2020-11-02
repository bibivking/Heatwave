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

case_name = (/"ERAI_SE39_Hyds_RUC","ERAI_SE39_Hyre"/)

path = "/g/data/w35/mm3972/model/wrf/NUWRF/LISWRF_configs/"

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; ____________________________ OPTIONS _______________________________
pic_name  = "ERAI_SE39_Hyre"

;______________________________ Data _________________________________
year_s = 2014
year_e = 2018
case_sum = dimsizes(case_name)
year_sum = year_e-year_s+1

if (year_s .eq. 1979) .and. (year_e .eq. 2018) then
  time_tot = (year_sum-1)*365 + 10
else
  time_tot = year_sum*365 + 2
end if

Qle      = new((/case_sum,time_tot/),float)
Qh       = new((/case_sum,time_tot/),float)
Qg       = new((/case_sum,time_tot/),float)
Evap     = new((/case_sum,time_tot/),float)
TVeg     = new((/case_sum,time_tot/),float)
ESoil    = new((/case_sum,time_tot/),float)
Qs       = new((/case_sum,time_tot/),float)
Qsb      = new((/case_sum,time_tot/),float)
Qrec     = new((/case_sum,time_tot/),float)

GWwb     = new((/case_sum,time_tot/),float)
Fwsoil   = new((/case_sum,time_tot/),float)
WatTab   = new((/case_sum,time_tot/),float)
SoilTemp = new((/case_sum,time_tot,6/),float)
SoilMoist= new((/case_sum,time_tot,6/),float)

Qle      = 0.
Qh       = 0.
Qg       = 0.
Evap     = 0.
TVeg     = 0.
ESoil    = 0.
Qs       = 0.
Qsb      = 0.
Qrec     = 0.
GWwb     = 0.
Fwsoil   = 0.
WatTab   = 0.
SoilTemp = 0.
SoilMoist= 0.

do case_num = 0,case_sum -1
  i = 0
  do year = year_s,year_e

      print("year = " +year)
      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if

      if year .eq. 1979 then
        month_s = 12
        month_e = 12
      else
        month_s = 1
        month_e = 12
      end if

      ; read LIS-CABLE
      do month = month_s, month_e
        print("month = "+month)
        if month .lt. 10 then
          filename = path+case_name(case_num)+"/LIS_output/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
        else
          filename = path+case_name(case_num)+"/LIS_output/LIS.CABLE."+year+month+"0100.d01.nc"
        end if
        f = addfile (filename,"r")
        ; ========================= READING =============================

        do day = 1, dom(month-1)
          Qle(case_num,i)      = avg(f->Qle_tavg(day-1,:,:))
          Qh(case_num,i)       = avg(f->Qh_tavg(day-1,:,:))
          Qg(case_num,i)       = avg(f->Qg_tavg(day-1,:,:))
          Evap(case_num,i)     = avg(f->Evap_tavg(day-1,:,:))*3600*24
          TVeg(case_num,i)     = avg(f->TVeg_tavg(day-1,:,:))*3600*24
          ESoil(case_num,i)    = avg(f->ESoil_tavg(day-1,:,:))*3600*24
          Qs(case_num,i)       = avg(f->Qs_tavg(day-1,:,:))*3600*24
          Qsb(case_num,i)      = avg(f->Qsb_tavg(day-1,:,:))*3600*24
          Qrec(case_num,i)     = avg(f->Qrec_tavg(day-1,:,:))*3600*24
          GWwb(case_num,i)     = avg(f->GWwb_inst(day-1,:,:))
          Fwsoil(case_num,i)   = avg(f->fwsoil_tavg(day-1,:,:))
          WatTab(case_num,i)   = avg(f->WaterTableD_tavg(day-1,:,:))/1000.

          do soil = 0,5
            SoilTemp(case_num,i,soil) = avg(f->SoilTemp_inst(day-1,soil,:,:))-273.15
            SoilMoist(case_num,i,soil)= avg(f->SoilMoist_inst(day-1,soil,:,:))
          end do
          i = i + 1
        end do ; day
      end do ; month
  end do ; year
  print("prepare to plot")
end do ; case_num

;_____________________________ PLOT ___________________________________
pic = "./plots/Spinup_LIS_"+pic_name+"_"+year_s+"-"+year_e

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@gsnMaximize        = True
res@vpHeightF          = 0.6
res@vpWidthF           = 0.6
res@vpXF               = 0.3
res@vpYF               = 0.85
res@pmLegendWidthF     = 0.14                     ; Resize legend width
res@pmLegendHeightF    = 0.11                     ; and height

res@pmLegendDisplayMode = "Always"
res@xyLineThicknesses  = (/  1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0/)          ; make second line thicker
res@xyDashPatterns     = (/ 0, 0, 0, 0, 0, 0, 0/)
res@xyLineColors       = (/"blue","red","orange","gray","black","green","pink"/)          ; change line color
res@xyExplicitLegendLabels = case_name
;res@trYMinF = 0


x             = ispan(1,time_tot,1)
res@tmXBMode  = "Explicit"
res@tmXBValues= ispan(31,time_tot,365)
res@tmXBLabels= tostring(ispan(year_s,year_e,1))

print(dimsizes(ispan(31,time_tot,365)))
print(dimsizes(ispan(year_s+1,year_e,1)))

; print("plot rainfall")
; ; Rain_rc    = regline(x, Rain)          ; slope
; ; Rain_arry  = new ( (/2,time_tot/), typeof(Rain))
; ; Rain_arry(0,:) = Rain
; ; Rain_arry(1,:) = Rain_rc*x + Rain_rc@yintercept
; if hydr_param then
;   res@tiYAxisString  = f->TotalPrecip_tavg@long_name+" (mm/d)"
; else
;   res@tiYAxisString  = f->Rainf_f_inst@long_name+" (mm/d)"
; end if
; res@tiXAxisString  = "days"
; plot = gsn_csm_xy(wks,x,Rain,res)
;
print("plot latent heat")
res@trYMaxF = 300.
res@trYMinF = -100.
; Qle_rc    = regline(x, Qle)          ; slope
; Qle_arry  = new ( (/2,time_tot/), typeof(Qle))
; Qle_arry(0,:) = Qle
; Qle_arry(1,:) = Qle_rc*x + Qle_rc@yintercept
res@tiYAxisString  = f->Qle_tavg@long_name+" ("+f->Qle_tavg@units+")"
res@tiXAxisString  = "days"
plot = gsn_csm_xy(wks,x,Qle,res)

print("plot sensible heat")
res@trYMaxF = 300.
res@trYMinF = -100.
; Qh_rc    = regline(x, Qh)          ; slope
; Qh_arry  = new ( (/2,time_tot/), typeof(Qh))
; Qh_arry(0,:) = Qh
; Qh_arry(1,:) = Qh_rc*x + Qh_rc@yintercept
res@tiYAxisString  = f->Qh_tavg@long_name+" ("+f->Qh_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qh,res)

print("plot ground heat")
res@trYMaxF = 100.
res@trYMinF = -100.
; Qg_rc    = regline(x, Qg)          ; slope
; Qg_arry  = new ( (/2,time_tot/), typeof(Qg))
; Qg_arry(0,:) = Qg
; Qg_arry(1,:) = Qg_rc*x + Qg_rc@yintercept
res@tiYAxisString  = f->Qg_tavg@long_name+" ("+f->Qg_tavg@units+")"
plot = gsn_csm_xy(wks,x,Qg,res)

print("plot ET")
res@trYMaxF = 8.
res@trYMinF = 0.
; Evap_rc    = regline(x, Evap)          ; slope
; Evap_arry  = new ( (/2,time_tot/), typeof(Evap))
; Evap_arry(0,:) = Evap
; Evap_arry(1,:) = Evap_rc*x + Evap_rc@yintercept
res@tiYAxisString  = f->Evap_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Evap,res)

print("plot trans")
res@trYMaxF = 5.
res@trYMinF = 0.
; TVeg_rc    = regline(x, TVeg)          ; slope
; TVeg_arry  = new ( (/2,time_tot/), typeof(TVeg))
; TVeg_arry(0,:) = TVeg
; TVeg_arry(1,:) = TVeg_rc*x + TVeg_rc@yintercept
res@tiYAxisString  = f->TVeg_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,TVeg,res)

print("plot esoil")
res@trYMaxF = 5.
res@trYMinF = 0.
; ESoil_rc    = regline(x, ESoil)          ; slope
; ESoil_arry  = new ( (/2,time_tot/), typeof(ESoil))
; ESoil_arry(0,:) = ESoil
; ESoil_arry(1,:) = ESoil_rc*x + ESoil_rc@yintercept
res@tiYAxisString  = f->ESoil_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,ESoil,res)

print("plot surface runoff")
res@trYMaxF = 5.
res@trYMinF = 0.
; Qs_rc    = regline(x, Qs)          ; slope
; Qs_arry  = new ( (/2,time_tot/), typeof(Qs))
; Qs_arry(0,:) = Qs
; Qs_arry(1,:) = Qs_rc*x + Qs_rc@yintercept
res@tiYAxisString  = f->Qs_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qs,res)

print("plot subsurface runoff")
res@trYMaxF = 5.
res@trYMinF = 0.
; Qsb_rc    = regline(x, Qsb)          ; slope
; Qsb_arry  = new ( (/2,time_tot/), typeof(Qsb))
; Qsb_arry(0,:) = Qsb
; Qsb_arry(1,:) = Qsb_rc*x + Qsb_rc@yintercept
res@tiYAxisString  = f->Qsb_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qsb,res)

print("plot aquifer recharge")
res@trYMaxF = 1.
res@trYMinF = -1.
; Qrec_rc    = regline(x, Qrec)          ; slope
; Qrec_arry  = new ( (/2,time_tot/), typeof(Qrec))
; Qrec_arry(0,:) = Qrec
; Qrec_arry(1,:) = Qrec_rc*x + Qrec_rc@yintercept
res@tiYAxisString  = f->Qrec_tavg@long_name+" (mm/d)"
plot = gsn_csm_xy(wks,x,Qrec,res)

print("plot fwsoil")
res@trYMaxF = 1.2
res@trYMinF = 0.
; fwsoil_rc    = regline(x, Fwsoil)          ; slope
; fwsoil_arry  = new ( (/2,time_tot/), typeof(Fwsoil))
; fwsoil_arry(0,:) = Fwsoil
; fwsoil_arry(1,:) = fwsoil_rc*x + fwsoil_rc@yintercept
res@tiYAxisString  = f->fwsoil_tavg@long_name+" ("+f->fwsoil_tavg@units+")"
plot = gsn_csm_xy(wks,x,Fwsoil,res)

print("plot water table")
res@trYMaxF = 50.
res@trYMinF = 0.
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

do soil = 0,5
  print("plot soil temp")
  res@trYMaxF = 35.
  res@trYMinF = 0.
  ; SoilTemp_rc    = regline(x, SoilTemp(:,soil))          ; slope
  ; SoilTemp_arry  = new ( (/2,time_tot/), typeof(SoilTemp))
  ; SoilTemp_arry(0,:) = SoilTemp(:,soil)
  ; SoilTemp_arry(1,:) = SoilTemp_rc*x + SoilTemp_rc@yintercept
  res@tiYAxisString  = "Layer "+(soil+1)+" Soil Temp (C)"
  plot = gsn_csm_xy(wks,x,SoilTemp(:,:,soil),res)
  ; delete(SoilTemp_rc)
  ; delete(SoilTemp_arry)
  delete(res@trYMaxF)
  delete(res@trYMinF)
end do
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
print("plot aquifer moisture")
res@trYMaxF = 0.5
res@trYMinF = 0.
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

end
