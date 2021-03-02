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

; ================================== Options ===================================
scale = "SE-AU"

if scale .eq. "SE-AU" then
  A = -40.
  B = -28.
  C = 140.
  D = 154.
end if

filename_rain = "/g/data/w35/Shared_data/Observations/AWAP_all_variables/monthly/rain/AWAP_monthly_rain_1970_2019.nc"
filename_tmax = "/g/data/w35/Shared_data/Observations/AWAP_all_variables/monthly/tmax/AWAP_monthly_tmax_1970_2019.nc"
filename_tmin = "/g/data/w35/Shared_data/Observations/AWAP_all_variables/monthly/tmin/AWAP_monthly_tmin_1970_2019.nc"

frain         = addfile(filename_rain,"r")
ftmax         = addfile(filename_tmax,"r")
ftmin         = addfile(filename_tmin,"r")

rain          = month_to_annual( dim_avg_n_Wrap(frain->rain(:,{A:B},{C:D}),(/1,2/)), 0 )
tmean         = month_to_annual(                                             \
                ( dim_avg_n_Wrap(ftmax->tmax(:,{A:B},{C:D}),(/1,2/))         \
                + dim_avg_n_Wrap(ftmin->tmin(:,{A:B},{C:D}),(/1,2/)) )/2., 1)
tmean_annual  = (dim_avg_n_Wrap(ftmax->tmax(:,{A:B},{C:D}),(/1,2/)) + dim_avg_n_Wrap(ftmin->tmin(:,{A:B},{C:D}),(/1,2/)) )/2.
tmean_annual!0 = "time"
tmean_sum     = month_to_seasonN(tmean_annual, (/"DJF"/))

print(rain)
print(tmean)
print(tmean_sum)
end
