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

path1 = "/g/data/w35/mm3972/model/cable/runs/test_mpi_serial/"

filename1 = path1+"cable_out_serial.nc"
filename2 = path1+"cable_out_mpi.nc"
filename3 = path1+"cable_out_serial_monthly.nc"
filename4 = path1+"cable_out_mpi_monthly.nc"

f1 = addfile (filename1,"r")
f2 = addfile (filename2,"r")
f3 = addfile (filename3,"r")
f4 = addfile (filename4,"r")

Fwsoil_all = new((/2,70128/),float)
Fwsoil_all(0,:) = f1->Fwsoil(:,0,0)
Fwsoil_all(1,:) = f2->Fwsoil(:,0,0)

printVarSummary(f3->Fwsoil)
printVarSummary(f4->Fwsoil)

Fwsoil_mon = new((/2,44/),float)
Fwsoil_mon(0,:) = f3->Fwsoil(0:43,0,0)
Fwsoil_mon(1,:) = f4->Fwsoil(0:43,0,0)

;_____________________________ PLOT ___________________________________

wks = gsn_open_wks("pdf","serial_vs_mpi")
gsn_define_colormap(wks,"BlueYellowRed")
res                    = True
res@gsnMaximize        = True
; res@vpHeightF          = 0.6
; res@vpWidthF           = 0.6
; res@vpXF               = 0.3
; res@vpYF               = 0.85
; res@pmLegendWidthF     = 0.14                     ; Resize legend width
; res@pmLegendHeightF    = 0.11                     ; and height

res@xyLineThicknesses  = (/  1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0/)          ; make second line thicker
res@pmLegendDisplayMode = "Always"
res@xyLineColors       = (/"blue","red","green","orange","red","red","red"/)          ; change line color
res@xyExplicitLegendLabels = (/"serial","mpi","GW_1-60","NoAQ"/)
;


x                  = ispan(1,70128,1)
res@trYMaxF = 2.5
res@trYMinF = 0.0
res@tiXAxisString  = "time step"
plot = gsn_csm_xy(wks,x,Fwsoil_all,res)

delete(x)
x                  = ispan(1,44,1)
res@tiXAxisString  = "month"
plot = gsn_csm_xy(wks,x,Fwsoil_mon,res)
end
