;----------------------------------------------------------------------
; wrf_gsn_5.ncl
;----------------------------------------------------------------------
; Concepts illustrated:
;   - Using gsn_csm scripts to plot WRF-ARW data
;   - Overlaying line contours, filled contours, and vectors on a map
;   - Setting the correct WRF map projection using wrf_map_resources
;   - Setting lots of resources to customize a WRF plot
;----------------------------------------------------------------------
; This script is meant to show the difference between plotting WRF
; data using wrf_xxx scripts, and using gsn_csm_xxx scripts.
;
; See wrf_nogsn_5.ncl for an example of using wrf_xxxx scripts to
; plot WRF data.
;----------------------------------------------------------------------
; In NCL Versions 6.3.1 and earlier, you will get these warnings which
; you can safely ignore:
;
; warning:start_lat is not a valid resource in wrf_gsn_contour at this time
; warning:start_lon is not a valid resource in wrf_gsn_contour at this time
; warning:end_lat is not a valid resource in wrf_gsn_contour at this time
; warning:end_lon is not a valid resource in wrf_gsn_contour at this time
; warning:mpNestTime is not a valid resource in map at this time
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;---Open WRF output file
  path = "/g/data/w35/mm3972/model/wrf/NUWRF/nuwrf_25km/LISWRF_configs/WRF_output/"

  ;filename = path+"wrfout_d01_2009-01-01_00:00:00" ;- Jan
  filename = path+"wrfout_d01_2009-02-01_03:00:00" ;- Feb
  ;filename = path+"wrfout_d01_2009-03-01_03:00:00"

  a        = addfile(filename,"r")
  time     = wrf_user_getvar(a,"times",-1) ; Times in file

  wks = gsn_open_wks("pdf","Feb-2009");"Feb-2009") ;"Jan-2009")

do it = 0, 223 ;247 -Mar ;223 -Feb  ; 248-Jan
  ;---Read several WRF variables at first time step
  ; it: Time in file to retrieve. -1 will retrieve all times in the file(s).
  slp = wrf_user_getvar(a,"slp",it)   ; sea level pressure
  wrf_smooth_2d( slp, 3 )             ; smooth slp
  tc  = wrf_user_getvar(a,"tc",it)    ; 3D temperature
  u   = wrf_user_getvar(a,"ua",it)    ; 3D U at mass points
  v   = wrf_user_getvar(a,"va",it)    ; 3D V at mass points

;---Now get the lowest (bottommost) level
  nl  = 0
  tc2 = tc(nl,:,:)
  u10 = u(nl,:,:)
  v10 = v(nl,:,:)
  ;tf2 = 1.8*tc2+32.                    ; Convert temperature to Fahrenheit
  u10 = u10*1.94386                    ; Convert wind into knots
  v10 = v10*1.94386

;---Change the metadata
  tc2@description = "Surface Temperature"
  tc2@units       = "degC"
  u10@units       = "kts"
  v10@units       = "kts"



;---Set common resources for all plots
  res                = True
  res@gsnFrame       = False
  res@gsnDraw        = False
  res@gsnLeftString  = ""
  res@gsnRightString = ""

;---Necessary for contours to be overlaid correctly on WRF projection
  res@tfDoNDCOverlay   = True          ; Tell NCL you are doing a native plot
; res@tfDoNDCOverlay   = "NDCViewport" ; can use this in NCL V6.5.0 or later

;---Temperature filled contour plot
  tc_res                      = res
  tc_res@cnFillOn             = True
  tc_res@cnLevelSelectionMode = "ExplicitLevels"
  tc_res@cnLevels             = ispan(0,50,5)
  tc_res@lbLabelFontHeightF   = 0.015
  tc_res@lbOrientation        = "Vertical"
  tc_res@pmLabelBarOrthogonalPosF = -0.005

  contour_tc = gsn_csm_contour(wks,tc2,tc_res)

;---SLP line contour plot
  levels      = ispan(900,1100,4)
  info_string = "Sea level pressure contours from 900 to 1100 by 4"

  slp_res                                = res
  slp_res@cnLineColor                    = "NavyBlue"
  slp_res@cnLevelSelectionMode           = "ExplicitLevels"
  slp_res@cnLevels                       = levels
  slp_res@cnLineLabelBackgroundColor     = -1    ; transparent
  slp_res@cnLineThicknessF               = 2.5
  slp_res@cnHighLabelsOn                 = True
  slp_res@cnLowLabelsOn                  = True
  slp_res@cnHighLabelBackgroundColor     = -1
  slp_res@cnLowLabelBackgroundColor      = -1
  slp_res@cnInfoLabelString              = info_string
  slp_res@cnInfoLabelFontColor           = "NavyBlue"
  slp_res@cnInfoLabelPerimOn             = False

  contour_psl = gsn_csm_contour(wks,slp,slp_res)

;---Wind vector plot
  vec_res                  = res
  vec_res@vcMinDistanceF   = 0.02
  vec_res@vcRefLengthF     = 0.02
  vec_res@vcMinFracLengthF = 0.2
  vec_res@vcGlyphStyle     = "WindBarb"
  vec_res@vcRefAnnoOn      = False

  vector = gsn_csm_vector(wks,u10,v10,vec_res)

;---Map plot
  map_res               = True
  map_res@gsnFrame      = False
  map_res@gsnDraw       = False
  map_res@tiMainString  = time(it)
  map_res@mpOutlineOn   = True
  map_res@mpCountyLineColor = "white"
  map_res@mpCountyLineThicknessF = 2.0
  map_res@gsnLeftString = tc2@description + " (" + tc@units + ")~C~" + \
                          slp@description + " (" + slp@units + ")~C~" + \
                          "Wind (" + u10@units + ")"
  map_res@gsnLeftStringFontHeightF = 0.01

;---Set map resources based on projection on WRF output file
  map_res = wrf_map_resources(a,map_res)

  map = gsn_csm_map(wks,map_res)

;---Overlay plots on map and draw.
  overlay(map,contour_tc)
  overlay(map,contour_psl)
  overlay(map,vector)

  draw(map)   ; This will draw all overlaid plots and the map
  frame(wks)

  u   = wrf_user_getvar(a,"ua",it)    ; 3D U at mass points
  v   = wrf_user_getvar(a,"va",it)    ; 3D V at mass points
  delete(slp)
  delete(tc)
  delete(u)
  delete(v)
  delete(nl)
  delete(tc2)
  delete(u10)
  delete(v10)
  delete(res)
  delete(tc_res)
  delete(slp_res)
  delete(vec_res)
  delete(map_res)
end do ; it
end
