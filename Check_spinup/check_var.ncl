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

; ================================ Pre-define ==================================
path = "/g/data/w35/mm3972/model/cable/runs/mmy_trunk_merged/"
case_name = (/"serial","mpi"/)
case_sum  = dimsizes(case_name)
filename  = new(case_sum,string)

var_2D_name = (/ "SWdown","LWdown","Tair","Rainf","Snowf","Qair","Wind","PSurf", \
                 "CO2air","Qmom","Qle","Qh","Qg","Qs","Qsb","Evap","ECanop","TVeg",\
                 "ESoil","HVeg","HSoil","RnetSoil","NEE","BaresoilT","SWE","SnowMelt",\
                 "SnowT","SnowDepth","SWnet","LWnet","Rnet","Albedo","RadT",\
                 "Tscrn","Qscrn","VegT","CanT","Fwsoil","CanopInt","LAI",\
                 "Ebal","Wbal","AutoResp","LeafResp","HeteroResp","GPP",\
                 "NPP","WatTable","GWMoist","SatFrac","Qrecharge","Area"/)

var_3D_name = (/ "SoilMoist","SoilMoistIce","SoilTemp","SMP"/)

; ================================ Dimsizes ====================================
do case_num = 0, case_sum-1
  filename(case_num) = path+"cable_out_"+case_name(case_num)+".nc"
end do

f        = addfile (filename(0),"r")
time_sum = dimsizes(f->time)
delete(f)

var_2D_sum = dimsizes(var_2D_name)
var_3D_sum = dimsizes(var_3D_name)

; ================================== PLOT ======================================
wks = gsn_open_wks("pdf","serial_vs_mpi_merged")
gsn_define_colormap(wks,"BlueYellowRed")
res                    = True
res@gsnMaximize        = True
; res@vpHeightF          = 0.6
; res@vpWidthF           = 0.6
; res@vpXF               = 0.3
; res@vpYF               = 0.85
; res@pmLegendWidthF     = 0.14                     ; Resize legend width
; res@pmLegendHeightF    = 0.11                     ; and height

res@pmLegendDisplayMode = "Always"
res@xyLineThicknesses  = (/  1.0, 1.0/)          ; make second line thicker
res@xyDashPatterns     = (/0, 0/)
res@xyLineColors       = (/"blue","red"/)          ; change line color
res@xyExplicitLegendLabels = case_name


x                  = ispan(1,time_sum,1)
; res@trYMaxF        = 2.5
; res@trYMinF        = 0.0
res@tiXAxisString  = "time step"

; =============================== Read and Plot ================================
do var_num = 0, var_2D_sum -1
  var_2D = new((/case_sum, time_sum/),float)
  do case_num = 0,case_sum -1
    f = addfile(filename(case_num),"r")
    var_2D(case_num,:) = f->$var_2D_name(var_num)$(:,0,0)
    delete(f)
  end do

  res@tiYAxisString  = var_2D_name(var_num)
  plot = gsn_csm_xy(wks,x,runave(var_2D(1,:),1000,0)-runave(var_2D(0,:),1000,0),res)
  delete(var_2D)
  delete(plot)
end do

do var_num = 0, var_3D_sum -1
  var_3D = new((/case_sum, time_sum/),float)
  do case_num = 0,case_sum -1
    f = addfile(filename(case_num),"r")
    var_3D(case_num,:) = dim_avg_n_Wrap(f->$var_3D_name(var_num)$(:,:,0,0),1)
    delete(f)
  end do

  res@tiYAxisString  = var_3D_name(var_num)
  plot = gsn_csm_xy(wks,x,runave(var_3D(1,:),1000,0)-runave(var_3D(0,:),1000,0),res)
  delete(var_3D)
  delete(plot)
end do

end
