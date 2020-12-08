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

fname_out  = "./restart_nc/restart_1999.nc"
finout     = addfile (fname_out, "w")     ; open output file

print(avg(finout->froot(0,:)))
print(avg(finout->froot(1,:)))
print(avg(finout->froot(2,:)))
print(avg(finout->froot(3,:)))
print(avg(finout->froot(4,:)))
print(avg(finout->froot(5,:)))

finout->froot(0,:)    = where(finout->iveg .ge. 1 .and. finout->iveg .le. 4, 0.01969, finout->froot(0,:))
finout->froot(1,:)    = where(finout->iveg .ge. 1 .and. finout->iveg .le. 4, 0.05008, finout->froot(1,:))
finout->froot(2,:)    = where(finout->iveg .ge. 1 .and. finout->iveg .le. 4, 0.12090, finout->froot(2,:))
finout->froot(3,:)    = where(finout->iveg .ge. 1 .and. finout->iveg .le. 4, 0.25017, finout->froot(3,:))
finout->froot(4,:)    = where(finout->iveg .ge. 1 .and. finout->iveg .le. 4, 0.34949, finout->froot(4,:))
finout->froot(5,:)    = where(finout->iveg .ge. 1 .and. finout->iveg .le. 4, 0.20967, finout->froot(5,:))


print(avg(finout->froot(0,:)))
print(avg(finout->froot(1,:)))
print(avg(finout->froot(2,:)))
print(avg(finout->froot(3,:)))
print(avg(finout->froot(4,:)))
print(avg(finout->froot(5,:)))


end
