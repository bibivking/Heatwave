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
;
; filo = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/CTL/restarts-spinup200yr/"
; system("/bin/cp " + filo)    ; remove if exists

fname_in   = "./restart_nc/restart_2010.nc"
fin        = addfile (fname_in, "r")     ; open output file

fname_out  = "./restart_nc/restart_2018.nc"
finout     = addfile (fname_out, "w")     ; open output file

print(avg(finout->wb))
print(avg(finout->wbice))
; print(avg(finout->GWwb))

finout->wb    = fin->wb*1.
finout->wbice = fin->wbice*1.
; finout->GWwb  = fin->GWwb*1.

print(avg(finout->wb))
print(avg(finout->wbice))
; print(avg(finout->GWwb))
; finout->tss = fin->tss ; Combined soil/snow temperature
; finout->tgg = fin->tgg ; Average layer soil temperature
; finout->tggsn = fin->tggsn ; Average layer snow temperature

end
