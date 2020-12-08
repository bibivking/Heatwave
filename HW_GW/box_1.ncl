;*********************************************
; box_1.ncl
;
; Concepts illustrated:
;   - Drawing box plots
;   - Explicitly setting the X tickmark labels in a box plot
;
;*********************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*********************************************
begin
;**********************************************
; Create some fake data
;**********************************************
  yval = new((/3,8/),"float",-999.)
  yval(0,0) = -3.
  yval(0,1) = -1.
  yval(0,2) = 1.5
  yval(0,3) = 4.2
  yval(0,4) = 6.
  yval(0,5) = -999.
  yval(0,6) = -999.
  yval(0,7) = -999.

  yval(1,0) = -1.
  yval(1,1) = 0.
  yval(1,2) = 1.
  yval(1,3) = 2.5
  yval(1,4) = 4.
  yval(1,5) = -999.
  yval(1,6) = -999.
  yval(1,7) = -999.

  yval(2,0) = -1.5
  yval(2,1) = 0.
  yval(2,2) = .75
  yval(2,3) = 2.
  yval(2,4) = 6.5
  yval(2,5) = -999.
  yval(2,6) = -999.
  yval(2,7) = -999.

  x = (/-3., -1., 1./)
;**********************************************
; create plot
;**********************************************
  wks = gsn_open_wks("png","box")               ; send graphics to PNG file

  res            = True                         ; plot mods desired
  res@tmXBLabels = (/"Control","-2Xna","2Xna"/) ; labels for each box
  res@tiMainString = "Default Box Plot"
;***********************************************
; the function boxplot will accept three different
; resource lists. In this default example, we set
; two of them to False.
;**********************************************
  plot = boxplot(wks,x,yval,False,res,False)
  draw(wks)                                    ; boxplot does not call these
  frame(wks)                                   ; for you

end
