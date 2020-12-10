;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./boxplot_mmy.ncl"
;*******************************************

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

; ================================== Options ===================================
year_s       = 2017
year_e       = 2017

scale        = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
heatwave     = "hw-90th" ; "hw-90th" ; "all-day"
case_name    = (/"CTL-30x3+30yr/outputs", "FREE_DRAIN/outputs-spinup30x3yr"/);,"PUMPx10","PUMPx100"/)

iveg_s = 1
iveg_e = 4

; ================================= Predefine ==================================
case_sum  = dimsizes(case_name)
year_sum = year_e-year_s+1

day_sum  = 0
do year = year_s,year_e
  if ( mod(year,4) .eq. 0) then
     day_sum = day_sum + 366
  else
     day_sum = day_sum + 365
  end if
end do

pic_name  = "HW_boxplot_"+heatwave+"_"+scale+"_CTL_vs_FREE_DRAIN_"+year_s+"-"+year_e

; =================================== Path =====================================
path = new(case_sum,"string")
do case_num = 0, case_sum -1
  path(case_num) = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"+case_name(case_num);+"/outputs"
end do

; ================================ Spitial Scale ===============================
if scale .eq. "Global" then
  A = -90.
  B = 90.
  C = 0.
  D = 360.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -28.
  C = 140.
  D = 154.
else if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
end if
end if
end if

; ================================== Start =====================================
filename = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
f        = addfile (filename,"r")
lat_sum  = dimsizes(f->Rainf(0,{A:B},0))
lon_sum  = dimsizes(f->Rainf(0,0,{C:D}))
def_value = f->Rainf@_FillValue

; PFT mask
PFT_mask  = where((f->iveg({A:B},{C:D}) .ge. iveg_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_e), 1, 0)
PFT_mask_multi = new((/day_sum,dimsizes(PFT_mask(:,0)),dimsizes(PFT_mask(0,:))/),float)
do z = 0, day_sum-1
  PFT_mask_multi(z,:,:) = PFT_mask
end do
delete(PFT_mask)
delete(f)

Rnet     = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
TVeg     = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
Fwsoil   = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
Tair     = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
deltaT   = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
SMsurf   = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
SM64cm   = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)
; SoilMoist= new((/case_sum,day_sum,6,lat_sum,lon_sum/),float,def_value)
WatTable = new((/case_sum,day_sum,lat_sum,lon_sum/),float,def_value)

Rnet     = 0.
TVeg     = 0.
Fwsoil   = 0.
Tair     = 0.
deltaT   = 0.
SMsurf   = 0.
SM64cm   = 0.
; SoilMoist= 0.
WatTable = 0.

; ___________ READ DATA ___________
day_step_s = 0
do year = year_s,year_e
  print("year = " +year)

    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    yr_tot_day = sum(dom)
    day_step_e = day_step_s + yr_tot_day - 1
    do case_num = 0,case_sum-1
      print("case_num = " +case_num)
      filename1 = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
      f1 = addfile (filename1,"r")
      TVeg(case_num,day_step_s:day_step_e,:,:)     = f1->TVeg(:,{A:B},{C:D})*3600.*24.
      Fwsoil(case_num,day_step_s:day_step_e,:,:)   = f1->Fwsoil(:,{A:B},{C:D})
      WatTable(case_num,day_step_s:day_step_e,:,:) = f1->WatTable(:,{A:B},{C:D})
      Tair(case_num,day_step_s:day_step_e,:,:)     = f1->Tair(:,{A:B},{C:D})-273.15
      deltaT(case_num,day_step_s:day_step_e,:,:)     = f1->VegT(:,{A:B},{C:D}) - f1->Tair(:,{A:B},{C:D})

      SMsurf(case_num,day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022 \
                                                     + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
                                                     + f1->SoilMoist(:,2,{A:B},{C:D})*(0.1 - 0.022 - 0.058)) \
                                                     / 0.1
      SM64cm(case_num,day_step_s:day_step_e,:,:)   = ( f1->SoilMoist(:,0,{A:B},{C:D})*0.022  \
                                                      + f1->SoilMoist(:,1,{A:B},{C:D})*0.058 \
                                                      + f1->SoilMoist(:,2,{A:B},{C:D})*0.154 \
                                                      + f1->SoilMoist(:,3,{A:B},{C:D})*0.409)\
                                                      / 0.643
      Rnet(case_num,day_step_s:day_step_e,:,:)     = f1->Rnet(:,{A:B},{C:D})
      ; SoilMoist(case_num,day_step_s:day_step_e,:,:,:)= f1->SoilMoist(:,:,{A:B},{C:D})
      delete(f1)
    end do ; case_num
    if heatwave .ne. "all-day" then
       day_step_s = day_step_e + 1
    end if
end do ; year

; ___________ MASK ___________
; mask HW
print("mask HW")
if heatwave .eq. "hw-90th" then
  filename_hw = "/g/data/w35/mm3972/scripts/ehfheatwaves/HW_Event_Indicator_1970-2019.nc" ; summer HW only
  fhw         = addfile(filename_hw,"r")

  HW = new((/day_sum,lat_sum,lon_sum/),float)
  day_start = 0
  do year = 1970, year_s - 1
    if ( mod(year,4) .eq. 0) then
       day_start = day_start + 366
    else
       day_start = day_start + 365
    end if
  end do
  HW = fhw->event(day_start:day_start+day_sum-1,:,:)*1.
  delete(filename_hw)
  delete(fhw)

  do case_num = 0,case_sum-1
    print(case_name(case_num))
    ; all PFT
    TVeg(case_num,:,:,:)     = where(HW .eq. 1, TVeg(case_num,:,:,:), def_value)
    Fwsoil(case_num,:,:,:)   = where(HW .eq. 1, Fwsoil(case_num,:,:,:), def_value)
    WatTable(case_num,:,:,:) = where(HW .eq. 1, WatTable(case_num,:,:,:), def_value)
    Tair(case_num,:,:,:)     = where(HW .eq. 1, Tair(case_num,:,:,:), def_value)
    deltaT(case_num,:,:,:)   = where(HW .eq. 1, deltaT(case_num,:,:,:), def_value)
    SMsurf(case_num,:,:,:)   = where(HW .eq. 1, SMsurf(case_num,:,:,:), def_value)
    SM64cm(case_num,:,:,:)   = where(HW .eq. 1, SM64cm(case_num,:,:,:), def_value)
    Rnet(case_num,:,:,:)     = where(HW .eq. 1, Rnet(case_num,:,:,:), def_value)
  end do
  delete(HW)
end if

; mask PFT
print("mask PFT")
TVeg_tree     = TVeg
Fwsoil_tree   = Fwsoil
WatTable_tree = WatTable
Tair_tree     = Tair
deltaT_tree   = deltaT
SMsurf_tree   = SMsurf
SM64cm_tree   = SM64cm
Rnet_tree     = Rnet

do case_num = 0,case_sum-1
  TVeg_tree(case_num,:,:,:)     = where(PFT_mask_multi .eq. 1, TVeg(case_num,:,:,:), def_value)
  Fwsoil_tree(case_num,:,:,:)   = where(PFT_mask_multi .eq. 1, Fwsoil(case_num,:,:,:), def_value)
  WatTable_tree(case_num,:,:,:) = where(PFT_mask_multi .eq. 1, WatTable(case_num,:,:,:), def_value)
  Tair_tree(case_num,:,:,:)     = where(PFT_mask_multi .eq. 1, Tair(case_num,:,:,:), def_value)
  deltaT_tree(case_num,:,:,:)   = where(PFT_mask_multi .eq. 1, deltaT(case_num,:,:,:), def_value)
  SMsurf_tree(case_num,:,:,:)   = where(PFT_mask_multi .eq. 1, SMsurf(case_num,:,:,:), def_value)
  SM64cm_tree(case_num,:,:,:)   = where(PFT_mask_multi .eq. 1, SM64cm(case_num,:,:,:), def_value)
  Rnet_tree(case_num,:,:,:)     = where(PFT_mask_multi .eq. 1, Rnet(case_num,:,:,:), def_value)
end do
delete(PFT_mask_multi)

; ___________ ND TO 1D ___________
; to 1d and remove missing values

; all PFTs
mask_no_mis      = ind( .not. ismissing(ndtooned(Fwsoil(0,:,:,:))))
length           = dimsizes(mask_no_mis)
print(length)
TVeg_no_mis      = new((/case_sum,length/), float, def_value)
Fwsoil_no_mis    = new((/case_sum,length/), float, def_value)
WatTable_no_mis  = new((/case_sum,length/), float, def_value)
Tair_no_mis      = new((/case_sum,length/), float, def_value)
deltaT_no_mis    = new((/case_sum,length/), float, def_value)
SMsurf_no_mis    = new((/case_sum,length/), float, def_value)
SM64cm_no_mis    = new((/case_sum,length/), float, def_value)
Rnet_no_mis      = new((/case_sum,length/), float, def_value)

; trees
mask_tree_no_mis = ind( .not. ismissing(ndtooned(Fwsoil_tree(0,:,:,:))))
length_tree      = dimsizes(mask_tree_no_mis)
print(length_tree)
print( sum(where (.not. ismissing(ndtooned(TVeg_tree(0,:,:,:))), 1, 0)))
TVeg_tree_no_mis      = new((/case_sum,length_tree/), float, def_value)
Fwsoil_tree_no_mis    = new((/case_sum,length_tree/), float, def_value)
WatTable_tree_no_mis  = new((/case_sum,length_tree/), float, def_value)
Tair_tree_no_mis      = new((/case_sum,length_tree/), float, def_value)
deltaT_tree_no_mis    = new((/case_sum,length_tree/), float, def_value)
SMsurf_tree_no_mis    = new((/case_sum,length_tree/), float, def_value)
SM64cm_tree_no_mis    = new((/case_sum,length_tree/), float, def_value)
Rnet_tree_no_mis      = new((/case_sum,length_tree/), float, def_value)

do case_num = 0,case_sum-1
  tmp                         = ndtooned(TVeg(case_num,:,:,:))
  TVeg_no_mis(case_num,:)     = tmp(mask_no_mis)
  tmp                         = ndtooned(Fwsoil(case_num,:,:,:))
  Fwsoil_no_mis(case_num,:)   = tmp(mask_no_mis)
  tmp                         = ndtooned(WatTable(case_num,:,:,:))
  WatTable_no_mis(case_num,:) = tmp(mask_no_mis)
  tmp                         = ndtooned(Tair(case_num,:,:,:))
  Tair_no_mis(case_num,:)     = tmp(mask_no_mis)
  tmp                         = ndtooned(deltaT(case_num,:,:,:))
  deltaT_no_mis(case_num,:)   = tmp(mask_no_mis)
  tmp                         = ndtooned(SMsurf(case_num,:,:,:))
  SMsurf_no_mis(case_num,:)   = tmp(mask_no_mis)
  tmp                         = ndtooned(SM64cm(case_num,:,:,:))
  SM64cm_no_mis(case_num,:)   = tmp(mask_no_mis)
  tmp                         = ndtooned(Rnet(case_num,:,:,:))
  Rnet_no_mis(case_num,:)     = tmp(mask_no_mis)
  delete(tmp)
  tmp                           = ndtooned(TVeg_tree(case_num,:,:,:))
  TVeg_tree_no_mis(case_num,:)  = tmp(mask_tree_no_mis)
  tmp                           = ndtooned(Fwsoil_tree(case_num,:,:,:))
  Fwsoil_tree_no_mis(case_num,:)= tmp(mask_tree_no_mis)
  tmp                           = ndtooned(WatTable_tree(case_num,:,:,:))
  WatTable_tree_no_mis(case_num,:) = tmp(mask_tree_no_mis)
  tmp                           = ndtooned(Tair_tree(case_num,:,:,:))
  Tair_tree_no_mis(case_num,:)  = tmp(mask_tree_no_mis)
  tmp                           = ndtooned(deltaT_tree(case_num,:,:,:))
  deltaT_tree_no_mis(case_num,:)= tmp(mask_tree_no_mis)
  tmp                           = ndtooned(SMsurf_tree(case_num,:,:,:))
  SMsurf_tree_no_mis(case_num,:)= tmp(mask_tree_no_mis)
  tmp                           = ndtooned(SM64cm_tree(case_num,:,:,:))
  SM64cm_tree_no_mis(case_num,:)= tmp(mask_tree_no_mis)
  tmp                           = ndtooned(Rnet_tree(case_num,:,:,:))
  Rnet_tree_no_mis(case_num,:)  = tmp(mask_tree_no_mis)
  delete(tmp)
end do
delete(TVeg)
delete(Fwsoil)
delete(WatTable)
delete(Tair)
delete(deltaT)
delete(SMsurf)
delete(SM64cm)
delete(Rnet)
delete(TVeg_tree)
delete(Fwsoil_tree)
delete(WatTable_tree)
delete(Tair_tree)
delete(deltaT_tree)
delete(SMsurf_tree)
delete(SM64cm_tree)
delete(Rnet_tree)

; ___________ Fwsoil Class ___________
; classify by Fwsoil levels
TVeg_cls      = new((/case_sum*10, length/), float, def_value)
Fwsoil_cls    = new((/case_sum*10, length/), float, def_value)
WatTable_cls  = new((/case_sum*10, length/), float, def_value)
Tair_cls      = new((/case_sum*10, length/), float, def_value)
deltaT_cls    = new((/case_sum*10, length/), float, def_value)
SMsurf_cls    = new((/case_sum*10, length/), float, def_value)
SM64cm_cls    = new((/case_sum*10, length/), float, def_value)
Rnet_cls      = new((/case_sum*10, length/), float, def_value)

TVeg_tree_cls      = new((/case_sum*10, length_tree/), float, def_value)
Fwsoil_tree_cls    = new((/case_sum*10, length_tree/), float, def_value)
WatTable_tree_cls  = new((/case_sum*10, length_tree/), float, def_value)
Tair_tree_cls      = new((/case_sum*10, length_tree/), float, def_value)
deltaT_tree_cls    = new((/case_sum*10, length_tree/), float, def_value)
SMsurf_tree_cls    = new((/case_sum*10, length_tree/), float, def_value)
SM64cm_tree_cls    = new((/case_sum*10, length_tree/), float, def_value)
Rnet_tree_cls      = new((/case_sum*10, length_tree/), float, def_value)

counter       = 0
do class = 0., 0.91 ,0.1
  do case_num = 0, case_sum -1
      print("counter = "+ counter)
      print(sum(where(Fwsoil_no_mis(case_num,:) .ge. class .and. Fwsoil_no_mis(case_num,:) .lt. class+0.1, \
      1, 0)))
      TVeg_cls(counter,:) = \
               where(Fwsoil_no_mis(case_num,:) .ge. class .and. Fwsoil_no_mis(case_num,:) .lt. class+0.1, \
               TVeg_no_mis(case_num,:), def_value)
      Fwsoil_cls(counter,:) = \
             where(Fwsoil_no_mis(case_num,:) .ge. class .and. Fwsoil_no_mis(case_num,:) .lt. class+0.1, \
             Fwsoil_no_mis(case_num,:), def_value)
      WatTable_cls(counter,:) = \
              where(Fwsoil_no_mis(case_num,:) .ge. class .and. Fwsoil_no_mis(case_num,:) .lt. class+0.1, \
              WatTable_no_mis(case_num,:), def_value)
      Tair_cls(counter,:) = \
               where(Fwsoil_no_mis(case_num,:) .ge. class .and. Fwsoil_no_mis(case_num,:) .lt. class+0.1, \
               Tair_cls(case_num,:), def_value)
      deltaT_cls(counter,:) = \
              where(Fwsoil_no_mis(case_num,:) .ge. class .and. Fwsoil_no_mis(case_num,:) .lt. class+0.1, \
              deltaT_no_mis(case_num,:), def_value)
      SMsurf_cls(counter,:) = \
               where(Fwsoil_no_mis(case_num,:) .ge. class .and. Fwsoil_no_mis(case_num,:) .lt. class+0.1, \
               SMsurf_no_mis(case_num,:), def_value)
      SM64cm_cls(counter,:) = \
              where(Fwsoil_no_mis(case_num,:) .ge. class .and. Fwsoil_no_mis(case_num,:) .lt. class+0.1, \
              SM64cm_no_mis(case_num,:), def_value)
      Rnet_cls(counter,:) = \
               where(Fwsoil_no_mis(case_num,:) .ge. class .and. Fwsoil_no_mis(case_num,:) .lt. class+0.1, \
               Rnet_no_mis(case_num,:), def_value)

      TVeg_tree_cls(counter,:) = \
              where(Fwsoil_tree_no_mis(case_num,:) .ge. class .and. Fwsoil_tree_no_mis(case_num,:) .lt. class+0.1, \
              TVeg_tree_no_mis(case_num,:), def_value)
      Fwsoil_tree_cls(counter,:) = \
            where(Fwsoil_tree_no_mis(case_num,:) .ge. class .and. Fwsoil_tree_no_mis(case_num,:) .lt. class+0.1, \
            Fwsoil_tree_no_mis(case_num,:), def_value)
      WatTable_tree_cls(counter,:) = \
             where(Fwsoil_tree_no_mis(case_num,:) .ge. class .and. Fwsoil_tree_no_mis(case_num,:) .lt. class+0.1, \
             WatTable_tree_no_mis(case_num,:), def_value)
      Tair_tree_cls(counter,:) = \
              where(Fwsoil_tree_no_mis(case_num,:) .ge. class .and. Fwsoil_tree_no_mis(case_num,:) .lt. class+0.1, \
              Tair_tree_cls(case_num,:), def_value)
      deltaT_tree_cls(counter,:) = \
             where(Fwsoil_tree_no_mis(case_num,:) .ge. class .and. Fwsoil_tree_no_mis(case_num,:) .lt. class+0.1, \
             deltaT_tree_no_mis(case_num,:), def_value)
      SMsurf_tree_cls(counter,:) = \
              where(Fwsoil_tree_no_mis(case_num,:) .ge. class .and. Fwsoil_tree_no_mis(case_num,:) .lt. class+0.1, \
              SMsurf_tree_no_mis(case_num,:), def_value)
      SM64cm_tree_cls(counter,:) = \
             where(Fwsoil_tree_no_mis(case_num,:) .ge. class .and. Fwsoil_tree_no_mis(case_num,:) .lt. class+0.1, \
             SM64cm_tree_no_mis(case_num,:), def_value)
      Rnet_tree_cls(counter,:) = \
              where(Fwsoil_tree_no_mis(case_num,:) .ge. class .and. Fwsoil_tree_no_mis(case_num,:) .lt. class+0.1, \
              Rnet_tree_no_mis(case_num,:), def_value)
      counter = counter + 1
  end do
end do

delete(TVeg_no_mis)
delete(Fwsoil_no_mis)
delete(WatTable_no_mis)
delete(Tair_no_mis)
delete(deltaT_no_mis)
delete(SMsurf_no_mis)
delete(SM64cm_no_mis)
delete(Rnet_no_mis)

delete(TVeg_tree_no_mis)
delete(Fwsoil_tree_no_mis)
delete(WatTable_tree_no_mis)
delete(Tair_tree_no_mis)
delete(deltaT_tree_no_mis)
delete(SMsurf_tree_no_mis)
delete(SM64cm_tree_no_mis)
delete(Rnet_tree_no_mis)


; ___________ Calc Stats ___________
; calculate and plot the median, minimum value, maximum value, and the 25th and 75th percentiles
TVeg_stat      = new((/case_sum*10, 5/), float, def_value)
Fwsoil_stat    = new((/case_sum*10, 5/), float, def_value)
WatTable_stat  = new((/case_sum*10, 5/), float, def_value)
Tair_stat      = new((/case_sum*10, 5/), float, def_value)
deltaT_stat    = new((/case_sum*10, 5/), float, def_value)
SMsurf_stat    = new((/case_sum*10, 5/), float, def_value)
SM64cm_stat    = new((/case_sum*10, 5/), float, def_value)
Rnet_stat      = new((/case_sum*10, 5/), float, def_value)

do counter = 0,case_sum*10-1

  tmp            = TVeg_cls(counter,:)
  tmp@_FillValue = def_value
  mask_no_mis_1  = ind( .not. ismissing(tmp))
  dimt           = dimsizes(mask_no_mis_1)
  print("dimt = "+dimt)
  x25            = round(.25*dimt,3)-1     ; -1 to account for NCL indexing starting
  x75            = round(.75*dimt,3)-1     ; at 0
  delete(tmp)

  tmp = TVeg_cls(counter,mask_no_mis_1)
  qsort(tmp)
  TVeg_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = Fwsoil_cls(counter,mask_no_mis_1)
  qsort(tmp)
  Fwsoil_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = WatTable_cls(counter,mask_no_mis_1)
  qsort(tmp)
  WatTable_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = Tair_cls(counter,mask_no_mis_1)
  qsort(tmp)
  Tair_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = deltaT_cls(counter,mask_no_mis_1)
  qsort(tmp)
  deltaT_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = SMsurf_cls(counter,mask_no_mis_1)
  qsort(tmp)
  SMsurf_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = SM64cm_cls(counter,mask_no_mis_1)
  qsort(tmp)
  SM64cm_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = Rnet_cls(counter,mask_no_mis_1)
  qsort(tmp)
  Rnet_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)
  delete(mask_no_mis_1)
  print(Rnet_stat(counter,:))
end do

delete(TVeg_cls)
delete(Fwsoil_cls)
delete(WatTable_cls)
delete(Tair_cls)
delete(deltaT_cls)
delete(SMsurf_cls)
delete(SM64cm_cls)
delete(Rnet_cls)

; trees
TVeg_tree_stat      = new((/case_sum*10, 5/), float, def_value)
Fwsoil_tree_stat    = new((/case_sum*10, 5/), float, def_value)
WatTable_tree_stat  = new((/case_sum*10, 5/), float, def_value)
Tair_tree_stat      = new((/case_sum*10, 5/), float, def_value)
deltaT_tree_stat    = new((/case_sum*10, 5/), float, def_value)
SMsurf_tree_stat    = new((/case_sum*10, 5/), float, def_value)
SM64cm_tree_stat    = new((/case_sum*10, 5/), float, def_value)
Rnet_tree_stat      = new((/case_sum*10, 5/), float, def_value)

do counter = 0,case_sum*10-1

  tmp            = TVeg_tree_cls(counter,:)
  tmp@_FillValue = def_value
  mask_no_mis_1  = ind( .not. ismissing(tmp))
  dimt           = dimsizes(mask_no_mis_1)
  print("dimt = "+dimt)
  x25            = round(.25*dimt,3)-1     ; -1 to account for NCL indexing starting
  x75            = round(.75*dimt,3)-1     ; at 0
  delete(tmp)

  tmp = TVeg_tree_cls(counter,mask_no_mis_1)
  qsort(tmp)
  TVeg_tree_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = Fwsoil_tree_cls(counter,mask_no_mis_1)
  qsort(tmp)
  Fwsoil_tree_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = WatTable_tree_cls(counter,mask_no_mis_1)
  qsort(tmp)
  WatTable_tree_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = Tair_tree_cls(counter,mask_no_mis_1)
  qsort(tmp)
  Tair_tree_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = deltaT_tree_cls(counter,mask_no_mis_1)
  qsort(tmp)
  deltaT_tree_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = SMsurf_tree_cls(counter,mask_no_mis_1)
  qsort(tmp)
  SMsurf_tree_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = SM64cm_tree_cls(counter,mask_no_mis_1)
  qsort(tmp)
  SM64cm_tree_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)

  tmp = Rnet_tree_cls(counter,mask_no_mis_1)
  qsort(tmp)
  Rnet_tree_stat(counter,:) = (/min(tmp),tmp(x25),dim_median(tmp),tmp(x75),max(tmp)/)
  delete(tmp)
  delete(mask_no_mis_1)
  print(Rnet_tree_stat(counter,:))
end do

delete(TVeg_tree_cls)
delete(Fwsoil_tree_cls)
delete(WatTable_tree_cls)
delete(Tair_tree_cls)
delete(deltaT_tree_cls)
delete(SMsurf_tree_cls)
delete(SM64cm_tree_cls)
delete(Rnet_tree_cls)

total_box       = case_sum*10*2+9
TVeg_stats      = new((/total_box, 5/), float, def_value)
Fwsoil_stats    = new((/total_box, 5/), float, def_value)
WatTable_stats  = new((/total_box, 5/), float, def_value)
Tair_stats      = new((/total_box, 5/), float, def_value)
deltaT_stats    = new((/total_box, 5/), float, def_value)
SMsurf_stats    = new((/total_box, 5/), float, def_value)
SM64cm_stats    = new((/total_box, 5/), float, def_value)
Rnet_stats      = new((/total_box, 5/), float, def_value)

do counter = 0,9
  TVeg_stats(counter*5  , :) = TVeg_stat(counter*2  , :)
  TVeg_stats(counter*5+1, :) = TVeg_stat(counter*2+1, :)
  TVeg_stats(counter*5+2, :) = TVeg_tree_stat(counter*2  , :)
  TVeg_stats(counter*5+3, :) = TVeg_tree_stat(counter*2+1, :)

  Fwsoil_stats(counter*5  , :) = Fwsoil_stat(counter*2  , :)
  Fwsoil_stats(counter*5+1, :) = Fwsoil_stat(counter*2+1, :)
  Fwsoil_stats(counter*5+2, :) = Fwsoil_tree_stat(counter*2  , :)
  Fwsoil_stats(counter*5+3, :) = Fwsoil_tree_stat(counter*2+1, :)

  WatTable_stats(counter*5  , :) = WatTable_stat(counter*2  , :)
  WatTable_stats(counter*5+1, :) = WatTable_stat(counter*2+1, :)
  WatTable_stats(counter*5+2, :) = WatTable_tree_stat(counter*2  , :)
  WatTable_stats(counter*5+3, :) = WatTable_tree_stat(counter*2+1, :)

  Tair_stats(counter*5  , :) = Tair_stat(counter*2  , :)
  Tair_stats(counter*5+1, :) = Tair_stat(counter*2+1, :)
  Tair_stats(counter*5+2, :) = Tair_tree_stat(counter*2  , :)
  Tair_stats(counter*5+3, :) = Tair_tree_stat(counter*2+1, :)

  deltaT_stats(counter*5  , :) = deltaT_stat(counter*2  , :)
  deltaT_stats(counter*5+1, :) = deltaT_stat(counter*2+1, :)
  deltaT_stats(counter*5+2, :) = deltaT_tree_stat(counter*2  , :)
  deltaT_stats(counter*5+3, :) = deltaT_tree_stat(counter*2+1, :)

  SMsurf_stats(counter*5  , :) = SMsurf_stat(counter*2  , :)
  SMsurf_stats(counter*5+1, :) = SMsurf_stat(counter*2+1, :)
  SMsurf_stats(counter*5+2, :) = SMsurf_tree_stat(counter*2  , :)
  SMsurf_stats(counter*5+3, :) = SMsurf_tree_stat(counter*2+1, :)

  SM64cm_stats(counter*5  , :) = SM64cm_stat(counter*2  , :)
  SM64cm_stats(counter*5+1, :) = SM64cm_stat(counter*2+1, :)
  SM64cm_stats(counter*5+2, :) = SM64cm_tree_stat(counter*2  , :)
  SM64cm_stats(counter*5+3, :) = SM64cm_tree_stat(counter*2+1, :)

  Rnet_stats(counter*5  , :) = Rnet_stat(counter*2  , :)
  Rnet_stats(counter*5+1, :) = Rnet_stat(counter*2+1, :)
  Rnet_stats(counter*5+2, :) = Rnet_tree_stat(counter*2  , :)
  Rnet_stats(counter*5+3, :) = Rnet_tree_stat(counter*2+1, :)
end do

x = ispan(1,total_box,1)

;**********************************************
; create plot
;**********************************************
wks = gsn_open_wks("png",pic_name)   ;"./plots/"+            ; send graphics to PNG file

res            = True                         ; plot mods desired
; res@gsnDraw    = True  ; Don't draw plots
; res@gsnFrame   = True  ; ����
; res@tmXBLabels = (/"0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5",\
;                    "0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1."/) ; labels for each box
res@tmXBLabels = (/"0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", \
                  "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1."/) ; labels for each box

res@tmXBValues   = ispan(1,total_box,5)
res@tiMainString = "Default Box Plot"

; resources for polylines that draws the boxes
llres                   = True
llres@gsLineThicknessF  = 2.5                ; line thickness

; resources that control color and width of boxes
opti           = True
opti@boxWidth  = 1.
opti@boxColors = (/"blue","red","green","orange", "blue","red","green","orange",\
                   "blue","red","green","orange", "blue","red","green","orange",\
                   "blue","red","green","orange", "blue","red","green","orange",\
                   "blue","red","green","orange", "blue","red","green","orange",\
                   "blue","red","green","orange", "blue","red","green","orange"/)

plot = boxplot_mmy(wks,x,deltaT_stats,opti,res,llres)
draw(wks)                                    ; boxplot does not call these
frame(wks)                                   ; for you

end
