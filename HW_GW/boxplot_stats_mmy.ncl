; ********************************************************
; Aim:
; Author: Mengyuan Mu
; ********************************************************

undef("boxplot_stats_mmy")
function boxplot_stats_mmy (var[*][*][*]:numeric,def_value:numeric)
;
;    var[*][*][*] : input data
;    def_value    : default values
;
begin

  var@_FillValue = def_value

  ; ___________ nd to 1d ___________
  tmp            = ndtooned(var)
  maskout_def    = ind( .not. ismissing(tmp))
  print(dimsizes(maskout_def))
  if dimsizes(maskout_def) .ge. 50 then
    var_masked     = tmp(maskout_def)
    var_masked@_FillValue = def_value

    ; ___________ Calc Stats ___________
    ; calculate and plot the median, minimum value, maximum value, and the 25th and 75th percentiles
    var_stat      = new(5, float, def_value)
    dimt          = dimsizes(var_masked)
    print("dimt = "+dimt)
    qsort(var_masked)
    x25            = round(.25*dimt,3)-1     ; -1 to account for NCL indexing starting
    x75            = round(.75*dimt,3)-1     ; at 0
    deltaQ         = var_masked(x75) - var_masked(x25)
    Minimum        = var_masked(x25) - 1.5*deltaQ ; the lowest data point excluding any outliers.
    Maximum        = var_masked(x75) + 1.5*deltaQ ; the largest data point excluding any outliers.

    var_stat = (/Minimum,var_masked(x25),dim_median(var_masked),var_masked(x75),Maximum/)
  else
    var_stat = (/-9999.,-9999.,-9999.,-9999.,-9999./)
  end if
  return(var_stat)
end
