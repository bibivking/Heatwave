
;MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
; Author  : Mengyuan Mu
; Version : 2.0 (22.06.2021)"
; Email   : mu.mengyuan815@gmail.com
;WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

;*******************************************
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

; ================================= Predefine ==================================
pic_name  = "P-E"
case_name = (/"CTL-30x3+30yr/outputs-daily","FREE_DRAIN/outputs-spinup30x3yr"/)
mode      = "plot";"all" ;"plot" ; "calc"

PFT       = "all" ; "tree" ;"nontree"
scale     = "SE-AU" ;"45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"
GLEAM_version = "3.5a";"3.3a"
case_sum  = dimsizes(case_name)

year_s   = 2000
year_e   = 2019
year_sum = year_e-year_s+1

time_tot = 0
do year = year_s,year_e
  if mod(year,4) .eq. 0 then
    time_tot = time_tot + 366
  else
    time_tot = time_tot + 365
  end if
end do

month_tot    = (year_e-year_s+1)*12

; 2001 -2009
total = 0
do year = 2001,2009
  if mod(year,4) .eq. 0 then
    total = total + 366
  else
    total = total + 365
  end if
end do
month_total = (2009-2001+1)*12+1

; set PFT
iveg_tree_s    = 1
iveg_tree_e    = 4
iveg_nontree_s = 5
iveg_nontree_e = 17

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

if (mode .eq. "calc" .or. mode .eq. "all") then
    ; _____________________________ Calc Bin Stats ________________________________
              ; 2001, 2017
    ts_s    = (/ 366, 6210/)
              ; 2009, 2019
    ts_e    = (/3652, 7304/)

                ; 2001, 2017
    mth_s    = (/   12, 204/)
                ; 2009, 2018
    mth_e    = (/  119, 227/)

    ; =================================== Path =====================================
    path = new(case_sum,"string")
    do case_num = 0, case_sum -1
      path(case_num) = "/g/data/w35/mm3972/model/cable/runs/AWAP_pumping/"+case_name(case_num)
    end do

    ; ================================== SET Dim ===================================
    filename  = path(0)+"/cable_out_"+year_s+"_SE_Aus.nc"
    f         = addfile (filename,"r")

    lat_size  = dimsizes(f->y)
    lon_size  = dimsizes(f->x)
    CABLE_lat = f->y
    CABLE_lon = f->x
    default_value = f->Rainf@_FillValue
    lat_sum = dimsizes(f->Rainf(0,{A:B},0))
    lon_sum = dimsizes(f->Rainf(0,0,{C:D}))

    Rain     = new((/case_sum, time_tot, lat_sum, lon_sum/), float, default_value)
    Rain_mth = new((/case_sum, month_tot, lat_sum, lon_sum/), float, default_value)
    Evap     = new((/case_sum, time_tot, lat_sum, lon_sum/), float, default_value)
    Rain     = 0.
    Rain_mth = 0.
    Evap     = 0.

    ; PFT mask
    PFT_mask = where((f->iveg({A:B},{C:D}) .ge. iveg_tree_s) .and. (f->iveg({A:B},{C:D}) .le. iveg_tree_e), 1, f->iveg({A:B},{C:D}))
    PFT_mask = where((PFT_mask .ge. iveg_nontree_s) .and. (PFT_mask .le. iveg_nontree_e), 2, PFT_mask)

    delete(f)

    ; ================================= Read Data ==================================
    do case_num = 0,case_sum-1
      day_step_s = 0
      count      = 0
      do year = year_s,year_e
        print("=====year = " +year)
        if ( mod(year,4) .eq. 0) then
           dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
        else
           dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
        end if
        filename1  = path(case_num)+"/cable_out_"+year+"_SE_Aus.nc"
        f1         = addfile (filename1,"r")
        yr_tot_day = sum(dom)
        day_step_e = day_step_s + yr_tot_day - 1

        Rain(case_num,day_step_s:day_step_e,:,:) = f1->Rainf(:,{A:B},{C:D})
        Evap(case_num,day_step_s:day_step_e,:,:) = f1->Evap(:,{A:B},{C:D})
        day_s = 0
        do month = 0,11
          day_e = day_s+dom(month)-1
          Rain_mth(case_num,count,:,:) = dim_sum_n_Wrap(f1->Rainf(day_s:day_e,{A:B},{C:D}),0)
          count = count + 1
          day_s = day_e + 1
        end do
        day_step_s = day_step_e + 1
        delete(f1)
      end do ; year
    end do ; case_num
    Rain = Rain*3600*24
    Evap = Evap*3600*24
    Rain_mth = Rain_mth*3600*24

    if PFT .eq. "tree" then
      do case_num = 0, case_sum - 1
        do day = 0,time_tot -1
          Rain(case_num,day,:,:)   = where(PFT_mask .eq. 1, Rain(case_num,day,:,:), default_value)
          Evap(case_num,day,:,:)   = where(PFT_mask .eq. 1, Evap(case_num,day,:,:), default_value)
        end do
        do month = 0,month_tot -1
          Rain_mth(case_num,month,:,:) = where(PFT_mask .eq. 1, Rain_mth(case_num,month,:,:), default_value)
        end do
      end do
    else if PFT .eq. "nontree" then
      do case_num = 0, case_sum - 1
        do day = 0,time_tot -1
          Rain(case_num,day,:,:)   = where(PFT_mask .eq. 2, Rain(case_num,day,:,:), default_value)
          Evap(case_num,day,:,:)   = where(PFT_mask .eq. 2, Evap(case_num,day,:,:), default_value)
        end do
        do month = 0,month_tot -1
          Rain_mth(case_num,month,:,:) = where(PFT_mask .eq. 2, Rain_mth(case_num,month,:,:), default_value)
        end do
      end do
    end if
    end if

    ; ================================ READ GLEAM ==================================
    ; use 2000-2019
    if GLEAM_version .eq. "3.3a" then
      path_GLEAM    = "/g/data/w35/Shared_data/Observations/Global_ET_products/GLEAM_v3_3/3_3a/daily/";
      gleam_ET      = path_GLEAM+year_s+"/E_"+year_s+"_GLEAM_v3.3a.nc"
    else if GLEAM_version .eq. "3.5a" then
      path_GLEAM    = "/g/data/ua8/GLEAM_v3-5/v3-5a/daily/";
      gleam_ET      = path_GLEAM+year_s+"/E_"+year_s+"_GLEAM_v3.5a.nc"
    end if
    end if

    G_ET          = addfile(gleam_ET,"r")
    GLEAM_def_val = doubletofloat(G_ET->E@_FillValue)
    tmp           = G_ET->E(:,{C:D},{A:B})
    GLEAM_lat     = tmp&lat
    GLEAM_lon     = tmp&lon
    delete(tmp)

    Evap_GLEAM    = new((/time_tot, lat_sum, lon_sum/), float, GLEAM_def_val)
    Evap_GLEAM    = GLEAM_def_val

    day_step_s = 0
    do year = year_s,year_e
      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if
      yr_tot_day = sum(dom)
      if GLEAM_version .eq. "3.3a" then
        gleam_ET   = path_GLEAM+year+"/E_"+year+"_GLEAM_v3.3a.nc"
      else if GLEAM_version .eq. "3.5a" then
        gleam_ET   = path_GLEAM+year+"/E_"+year+"_GLEAM_v3.5a.nc"
      end if
      end if
      G_ET       = addfile (gleam_ET,"r")
      day_step_e = day_step_s + yr_tot_day - 1
      E_temp     = G_ET->E(:,{C:D}, {A:B})
      E_tmp      = doubletofloat(E_temp(time|:,lat|:, lon|:))
      delete(E_temp)

      ; Regrid
      Evap_GLEAM(day_step_s:day_step_e,:,:) = \
      linint2_Wrap(GLEAM_lon,GLEAM_lat,E_tmp,False,CABLE_lon,CABLE_lat,0)

      delete(gleam_ET)
      delete(G_ET)
      delete(E_tmp)
      day_step_s = day_step_e + 1
    end do
    Evap_GLEAM = Evap_GLEAM

    ; ___________________________ Mask GLEAM ______________________________
    if PFT .eq. "tree" then
      do day = 0,time_tot -1
        Evap_GLEAM(day,:,:) = where(PFT_mask .eq. 1 .and. .not. ismissing(Evap_GLEAM(day,:,:)),\
                                    Evap_GLEAM(day,:,:), default_value)
      end do
    else if PFT .eq. "nontree" then
      do day = 0,time_tot -1
        Evap_GLEAM(day,:,:) = where(PFT_mask .eq. 2 .and. .not. ismissing(Evap_GLEAM(day,:,:)),\
                                    Evap_GLEAM(day,:,:), default_value)
      end do
    end if
    end if
    Evap_GLEAM@_FillValue = default_value


    ; ================================ READ DOLCE ==================================
    ; use 2000-2019

    path_DOLCE    = "/g/data/w35/mm3972/data/DOLCE/v3";

    DOLCE_ET      = path_DOLCE+"/DOLCE_v3_"+year_s+".nc"
    D_ET          = addfile(DOLCE_ET,"r")
    lat_sum_DOLCE = dimsizes(D_ET->hfls(0,{A:B},0))
    lon_sum_DOLCE = dimsizes(D_ET->hfls(0,0,{C:D}))
    DOLCE_def_val = D_ET->hfls@_FillValue

    ; Values
    DOLCE         = new((/month_tot, lat_sum_DOLCE, lon_sum_DOLCE/), float, DOLCE_def_val)
    DOLCE         = 0.
    DOLCE!0       = "month"
    DOLCE!1       = "lat"
    DOLCE!2       = "lon"
    DOLCE&lat     = doubletofloat(D_ET->lat({A:B}))
    DOLCE&lon     = doubletofloat(D_ET->lon({C:D}))

    ; SD
    DOLCE_SD      = new((/month_tot, lat_sum_DOLCE, lon_sum_DOLCE/), float, DOLCE_def_val)
    DOLCE_SD      = 0.
    DOLCE_SD!0    = "month"
    DOLCE_SD!1    = "lat"
    DOLCE_SD!2    = "lon"
    DOLCE_SD&lat  = doubletofloat(D_ET->lat({A:B}))
    DOLCE_SD&lon  = doubletofloat(D_ET->lon({C:D}))

    count = 0

    do year = year_s,year_e-1
      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if
      DOLCE_ET    = path_DOLCE+"/DOLCE_v3_"+year+".nc"
      D_ET        = addfile(DOLCE_ET,"r")
      do month = 0,11
        DOLCE(count,:,:)     = D_ET->hfls(month,{A:B},{C:D})*dom(month)/28.94
        DOLCE_SD(count,:,:)  = D_ET->hfls_sd(month,{A:B},{C:D})*dom(month)/28.94
        count = count + 1
      end do
      delete(D_ET)
    end do
    print("=====READ DOLCE")
    printVarSummary(DOLCE)
    printVarSummary(DOLCE_SD)

    if (any(isnan_ieee(DOLCE))) then
      replace_ieeenan (DOLCE, default_value, 0)
    end if

    if (any(isnan_ieee(DOLCE_SD))) then
      replace_ieeenan (DOLCE_SD, default_value, 0)
    end if
    DOLCE@_FillValue    = default_value
    DOLCE_SD@_FillValue = default_value


    print("=====Change DOLCE default value to CABLE default")
    printVarSummary(DOLCE)
    printVarSummary(DOLCE_SD)

    ; Regrid DOLCE
    Evap_DOLCE    = linint2_Wrap(DOLCE&lon,DOLCE&lat,DOLCE,False,CABLE_lon,CABLE_lat,0)
    Evap_DOLCE_SD = linint2_Wrap(DOLCE_SD&lon,DOLCE_SD&lat,DOLCE_SD,False,CABLE_lon,CABLE_lat,0)
    delete(DOLCE)
    delete(DOLCE_SD)

    Evap_DOLCE!0       = "month"
    Evap_DOLCE!1       = "lat"
    Evap_DOLCE!2       = "lon"
    Evap_DOLCE&lat     = CABLE_lat
    Evap_DOLCE&lon     = CABLE_lon
    Evap_DOLCE@_FillValue = default_value

    Evap_DOLCE_SD!0    = "month"
    Evap_DOLCE_SD!1    = "lat"
    Evap_DOLCE_SD!2    = "lon"
    Evap_DOLCE_SD&lat  = CABLE_lat
    Evap_DOLCE_SD&lon  = CABLE_lon
    Evap_DOLCE_SD@_FillValue = default_value

    print("=====Check the default value of DOLCE after regrid")
    printVarSummary(Evap_DOLCE)
    printVarSummary(Evap_DOLCE_SD)

    ; ___________________________ Mask DOLCE ______________________________
    if PFT .eq. "tree" then
      do month = 0,month_tot -1
        Evap_DOLCE(month,:,:)    = where(PFT_mask .eq. 1 .and. .not. ismissing(Evap_DOLCE(month,:,:)),\
                                    Evap_DOLCE(month,:,:), default_value)
        Evap_DOLCE_SD(month,:,:) = where(PFT_mask .eq. 1 .and. .not. ismissing(Evap_DOLCE_SD(month,:,:)),\
                                    Evap_DOLCE_SD(month,:,:), default_value)
      end do
    else if PFT .eq. "nontree" then
      do month = 0,month_tot -1
        Evap_DOLCE(month,:,:)    = where(PFT_mask .eq. 2 .and. .not. ismissing(Evap_DOLCE(month,:,:)),\
                                    Evap_DOLCE(month,:,:), default_value)
        Evap_DOLCE_SD(month,:,:) = where(PFT_mask .eq. 2 .and. .not. ismissing(Evap_DOLCE_SD(month,:,:)),\
                                    Evap_DOLCE_SD(month,:,:), default_value)
      end do
    end if
    end if

    if (any(isnan_ieee(Evap_DOLCE))) then
      print("=====Error Evap_DOLCE has NaN")
    end if

    if (any(isnan_ieee(Evap_DOLCE_SD))) then
      print("=====Error Evap_DOLCE_SD has NaN")
    end if

    print("=====Check the default value of DOLCE after mask")
    printVarSummary(Evap_DOLCE)
    printVarSummary(Evap_DOLCE_SD)

    ; ___________________________Change Def Value__________________________________
    print("=====Change default value to -9999.")

    default_new_value = -9999.

    Evap_DOLCE = where(ismissing(Evap_DOLCE), default_new_value, Evap_DOLCE)
    Evap_DOLCE@_FillValue = default_new_value
    printVarSummary(Evap_DOLCE)

    Evap_DOLCE_SD = where(ismissing(Evap_DOLCE_SD), default_new_value, Evap_DOLCE_SD)
    Evap_DOLCE_SD@_FillValue = default_new_value
    printVarSummary(Evap_DOLCE_SD)

    Evap_GLEAM = where(ismissing(Evap_GLEAM), default_new_value, Evap_GLEAM)
    Evap_GLEAM@_FillValue = default_new_value
    printVarSummary(Evap_GLEAM)

    Rain = where(ismissing(Rain), default_new_value, Rain)
    Rain@_FillValue = default_new_value
    printVarSummary(Rain)

    Evap = where(ismissing(Evap), default_new_value, Evap)
    Evap@_FillValue = default_new_value
    printVarSummary(Evap)

    Rain_mth = where(ismissing(Rain_mth), default_new_value, Rain_mth)
    Rain_mth@_FillValue = default_new_value
    printVarSummary(Rain_mth)

    ; ============================= plot time series ===============================

    ; _________________ CABLE & GLEAM ___________________
    P         = dim_avg_n_Wrap(Rain,(/2,3/))
    E         = dim_avg_n_Wrap(Evap,(/2,3/))
    GLEAM_E   = dim_avg_n_Wrap(Evap_GLEAM,(/1,2/))
    P_E       = dim_avg_n_Wrap(Rain,(/2,3/))
    P_E       = P_E - dim_avg_n_Wrap(Evap,(/2,3/))
    P_GLEAM_E = dim_avg_n_Wrap(Rain(0,:,:,:),(/1,2/))
    P_GLEAM_E = P_GLEAM_E - dim_avg_n_Wrap(Evap_GLEAM,(/1,2/))

    P_E_acl   = new((/(case_sum+1)*2,total/),float, default_new_value)
    P_E_acl   = default_new_value

    ; _________________ DOLCE ___________________
    DOLCE_E   = dim_avg_n_Wrap(Evap_DOLCE,(/1,2/))
    P_MTH     = dim_avg_n_Wrap(Rain_mth,(/2,3/))
    P_DOLCE_E = dim_avg_n_Wrap(Rain_mth(0,:,:,:),(/1,2/))
    P_DOLCE_E = P_DOLCE_E - dim_avg_n_Wrap(Evap_DOLCE,(/1,2/))

    P_DOLCE_E_plus_SD   = dim_avg_n_Wrap(Rain_mth(0,:,:,:),(/1,2/))
    P_DOLCE_E_plus_SD   = P_DOLCE_E_plus_SD - dim_avg_n_Wrap(Evap_DOLCE,(/1,2/)) + dim_avg_n_Wrap(Evap_DOLCE_SD,(/1,2/))
    P_DOLCE_E_minus_SD  = dim_avg_n_Wrap(Rain_mth(0,:,:,:),(/1,2/))
    P_DOLCE_E_minus_SD  = P_DOLCE_E_minus_SD - dim_avg_n_Wrap(Evap_DOLCE,(/1,2/)) - dim_avg_n_Wrap(Evap_DOLCE_SD,(/1,2/))

    P_E_DOLCE_acl = new((/2*3,month_total/),float, default_new_value)
    P_E_DOLCE_acl = default_new_value

    P_E_DOLCE_acl_SD = new((/2*2,month_total/),float, default_new_value)
    P_E_DOLCE_acl_SD = default_new_value

    count = 0
    do cnt = 0, dimsizes(ts_s) -1
      do case_num = 0, case_sum -1
        P_E_acl(count,ts_s(cnt)-ts_s(cnt):ts_e(cnt)-ts_s(cnt)) = dim_cumsum_n_Wrap(P_E(case_num,ts_s(cnt):ts_e(cnt)), 0, 0)
        count = count + 1
      end do
      P_E_acl(count,ts_s(cnt)-ts_s(cnt):ts_e(cnt)-ts_s(cnt))   = dim_cumsum_n_Wrap(P_GLEAM_E(ts_s(cnt):ts_e(cnt)), 0, 0)
      count = count + 1
      P_E_DOLCE_acl(cnt*3,   mth_s(cnt)-mth_s(cnt)+1:mth_e(cnt)-mth_s(cnt)+1) = \
                                                               dim_cumsum_n_Wrap(P_DOLCE_E(mth_s(cnt):mth_e(cnt)), 0, 0)
      P_E_DOLCE_acl(cnt*3+1, mth_s(cnt)-mth_s(cnt)+1:mth_e(cnt)-mth_s(cnt)+1) = \
                                                               dim_cumsum_n_Wrap(P_DOLCE_E_plus_SD(mth_s(cnt):mth_e(cnt)), 0, 0)
      P_E_DOLCE_acl(cnt*3+2, mth_s(cnt)-mth_s(cnt)+1:mth_e(cnt)-mth_s(cnt)+1) = \
                                                               dim_cumsum_n_Wrap(P_DOLCE_E_minus_SD(mth_s(cnt):mth_e(cnt)), 0, 0)
    end do

    P_E_DOLCE_acl(0,0)    = 0.
    P_E_DOLCE_acl(3,0)    = 0.
    P_E_DOLCE_acl(1:2,0)  = 0.
    P_E_DOLCE_acl(4:5,0)  = 0.

    r_CTL_MD = escorc(P(0,ts_s(0):ts_e(0)), E(0,ts_s(0):ts_e(0)))
    r_CTL_RD = escorc(P(0,ts_s(1):ts_e(1)), E(0,ts_s(1):ts_e(1)))

    r_FD_MD = escorc(P(1,ts_s(0):ts_e(0)), E(1,ts_s(0):ts_e(0)))
    r_FD_RD = escorc(P(1,ts_s(1):ts_e(1)), E(1,ts_s(1):ts_e(1)))

    r_GLEAM_MD = escorc(P(0,ts_s(0):ts_e(0)), GLEAM_E(ts_s(0):ts_e(0)))
    r_GLEAM_RD = escorc(P(0,ts_s(1):ts_e(1)), GLEAM_E(ts_s(1):ts_e(1)))

    print("=====r_CTL_MD: "+r_CTL_MD+"  r_CTL_RD: "+r_CTL_RD)
    print("=====r_FD_MD: "+r_FD_MD+"  r_FD_RD: "+r_FD_RD)
    print("=====r_GLEAM_MD: "+r_GLEAM_MD+"  r_GLEAM_RD: "+r_GLEAM_RD)

    ; __________________________ WRITE TXT OUTPUT _____________________________
    print("=====txt output")
    printVarSummary(P_E_acl)
    printVarSummary(P_E_DOLCE_acl)

    option        = True
    option@row    = False

    option@fout   =  "./txt/Fig1b_time_accumulate_P-E_GLEAM_"+GLEAM_version+"_PFT-"+PFT+"_"+year_s+"-"+year_e+".txt"
    format        = total+"f20.4"
    write_matrix (P_E_acl, format, option)
    delete(option@fout)


    option@fout   =  "./txt/Fig1b_time_accumulate_P-E_DOLCE_PFT-"+PFT+"_"+year_s+"-"+year_e+".txt"
    format        = month_total+"f20.4"
    write_matrix (P_E_DOLCE_acl, format, option)
    delete(option@fout)


else if (mode .eq. "plot" .or. mode .eq. "all") then
    if mode .eq. "plot" then
      default_new_value = -9999.
      txt_name  =  "./txt/Fig1b_time_accumulate_P-E_GLEAM_"+GLEAM_version+"_PFT-"+PFT+"_"+year_s+"-"+year_e+".txt"
      P_E_acl   = asciiread(txt_name,(/6,total/),"float")
      P_E_acl@_FillValue = default_new_value
      txt_name        =  "./txt/Fig1b_time_accumulate_P-E_DOLCE_PFT-"+PFT+"_"+year_s+"-"+year_e+".txt"
      P_E_DOLCE_acl   = asciiread(txt_name,(/6,month_total/),"float")
      P_E_DOLCE_acl@_FillValue = default_new_value
      print(P_E_DOLCE_acl(0::3,:))
    end if

    ; ______________________________ Plotting ______________________________
    pic = "./plots/Fig1b_time_accumulate_P-E_GLEAM_"+GLEAM_version+"_PFT-"+PFT+"_"+year_s+"-"+year_e
    wks = gsn_open_wks("pdf",pic)
    gsn_define_colormap(wks,"rainbow") ;"ViBlGrWhYeOrRe") ;"BlueYellowRed")

    res                         = True
    res@gsnDraw                 = False             ; don't draw yet
    res@gsnFrame                = False             ; don't advance frame yet
    res@gsnMaximize             = False

    res@vpHeightF 	            = 0.4               ; change aspect ratio of plot
    res@vpWidthF 	              = 0.8
    ; res@vpXF                    = 0.1
    ; res@vpYF                    = 0.9
    res@txFontHeightF           = 0.015            ; change title font heights
    res@tiYAxisFontThicknessF   = 2.
    res@tmXBMajorThicknessF     = 2.
    res@tmXBMinorThicknessF     = 2.
    res@tmYLMajorThicknessF     = 2.
    res@tmYLMinorThicknessF     = 2.
    res@tmYLLabelFontThicknessF = 2.
    res@tmBorderThicknessF      = 2.
    ; res@txFont                  = 25

    res@lgPerimOn              = False      ; No legend perimeter.
    res@pmLegendParallelPosF    = 0.36       ; Move to right
    res@pmLegendOrthogonalPosF  = -1.27     ; Move into plot
    res@pmLegendWidthF          = 0.24                     ; resize legend width
    res@pmLegendHeightF         = 0.21                     ; and height
    res@lgBoxMinorExtentF       = 0.28       ; Shorten the legend lines

    ; res@tiXAxisFontHeightF      = 0.015
    ; res@tiYAxisFontHeightF      = 0.015
    ; res@tmXBLabelFontHeightF    = 0.015
    ; res@tmYLLabelFontHeightF    = 0.015
    ; res@gsnPaperOrientation     = "portrait";"portrait";"landscape" ;"auto"
    res@pmLegendDisplayMode     = "Never";"Always"

    res@trXMinF                 = 0
    res@trXMaxF                 = total
    res@trYMinF                 = -110.
    res@trYMaxF                 = 800.

    res@tmXBMode                = "Explicit"
    res@tmXBValues              = ispan(1,total,366)
    res@tmXBLabels              = (/"2001~C~2017","","2003~C~2019","","2005","","2007","","2009"/)
    x                           = ispan(1,total,1)

    res@xyLineThicknesses      = (/ 3.5, 3.5, 3.5,  3.5, 3.5, 3.5, 3.5, 3.5/)
    res@xyDashPatterns         = (/  0,  0,  0,   0,  0,  0,  0,  0/)
    res@xyLineOpacities        = (/ 1., 1., 1.,  1., 1., 1., 1., 1./) ;0.5, 0.5, 0.5/)
    ; res@xyLineColors           = (/"darkgreen","orange", "black","blue3",  "green3","gold","gray50", "deepskyblue"  /)
    ; res@xyExplicitLegendLabels = (/"GW MD (r="+sprintf("%4.2f",r_CTL_MD)+")", \
    ;                                "FD MD  (r="+sprintf("%4.2f",r_FD_MD)+")",\
    ;                                "GLEAM MD", \
    ;                                "DOLCE MD", \
    ;                                "GW RD (r="+sprintf("%4.2f",r_CTL_RD)+")",\
    ;                                "FD RD  (r="+sprintf("%4.2f",r_FD_RD)+")",\
    ;                                "GLEAM RD", \
    ;                                "DOLCE RD"  /)  ; create explicit labels
    ; res@lgItemOrder            = (/5, 7, 3, 1, 4, 6, 2, 0/)
    ; res@lgJustification        = "TopRight"


    ; res@tiXAxisString        = "months"
    ;**************** ͼ�������� *****************
    pres                            = True
    pres@wkWidth                    = 2500
    pres@wkHeight                   = 2500
    pres@gsnMaximize                = False
    pres@gsnPanelYWhiteSpacePercent = 0.1
    pres@gsnPanelXWhiteSpacePercent = 0.1
    ; pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
    ; pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
    ; pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
    ; pres@lbTitleFont        = 25
    ; pres@lbLabelFont        = 25
    ; pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
    ; pres@lbTitleFontHeightF = 0.2
    ; pres@lbLabelFontHeightF = 0.2

    ; -------------- Fig 1b GW, FD, GLEAM lines --------------
    print("=====plot GW,FD,GLEAM")
    print("=====dimsizes(x) ="+dimsizes(x))
    res@tiYAxisString = "Accumulated P-E (mm)"

    res@xyLineColors  = (/"darkgreen","orange", "black", "green3","gold","gray50" /)
    plot              = gsn_csm_xy(wks,x,P_E_acl,res)


    ; -------------- Fig 1b DOLCE lines  --------------
    print("=====plot DOLCE")
    res1              = True
    res1@trYMinF      = -110.
    res1@trYMaxF      = 800.
    res1@trXMinF      = 0
    res1@trXMaxF      = total

    y            =  new(109,integer)
    count        = 0
    y(0)         = 0

    do year = 2001,2009
      print("=====year = " +year)
      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if
      do mth = 0,11
        count    = count + 1
        y(count) = y(count-1)+dom(mth)
      end do
    end do

    print("=====dimsizes(y) ="+dimsizes(y))
    print(y)
    res1@xyLineThicknesses  = (/ 3.5, 3.5, 3.5,  3.5, 3.5, 3.5, 3.5, 3.5/)
    res1@xyDashPatterns     = (/  0,  0,  0,   0,  0,  0,  0,  0/)
    res1@xyLineOpacities    = (/ 1., 1., 1.,  1., 1., 1., 1., 1./) ;0.5, 0.5, 0.5/)
    res1@xyLineColors       = (/"blue3","deepskyblue" /)
    plot_DOLCE   = gsn_csm_xy(wks,y,P_E_DOLCE_acl(0::3,:),res1)
    overlay(plot,plot_DOLCE)
    delete(res1@xyLineThicknesses)
    delete(res1@xyDashPatterns)
    delete(res1@xyLineOpacities)
    delete(res1@xyLineColors)

    ; ; ________________ add error bar ________________
    ; ; gsresy = True
    ; ; gsresy@gsLineColor      = "blue"                       ; add y error bars
    ; ; gsresy@gsLineThicknessF = 2.0
    ;
    ; print("=====add uncertainty")
    ; res1@xyLineColor        = -1
    ; res1@gsnXYFillOpacities = 0.1
    ;
    ; res1@gsnXYFillColors    = "blue3"
    ; plot_uncert             = gsn_csm_xy (wks,y,P_E_DOLCE_acl(1:2,:),res1)
    ; overlay(plot,plot_uncert)
    ; res1@gsnXYFillColors    = "deepskyblue"
    ; plot_uncert             = gsn_csm_xy (wks,y,P_E_DOLCE_acl(4:5,:),res1)
    ; overlay(plot,plot_uncert)


    ; ________________ add order name ________________
    txres                  = True
    txres@txPerimOn        = False
    txres@txFont           = 25
    txres@txFontHeightF    = 0.025
    amres                  = True
    amres@amParallelPosF   = -0.487    ; This is the right edge of the plot.
    amres@amOrthogonalPosF = -0.46    ; This is the bottom edge of the plot.
    amres@amJust           = "TopLeft"

    txid1   = gsn_create_text(wks, "(b)", txres)
    annoid1 = gsn_add_annotation(plot, txid1, amres)

    ; ________________ add legend ________________
    genres                           = True
    genres@XPosPercent               = 10                      ; orientation on page
    genres@YPosPercent               = 92
    genres@ItemSpacePercent          = 6
    textres                          = True
    ; textres@lgLabels                 = (/"GW MD (r="+sprintf("%4.2f",r_CTL_MD)+")", \
    ;                                      "FD MD  (r="+sprintf("%4.2f",r_FD_MD)+")",\
    ;                                      "GLEAM MD", \
    ;                                      "DOLCE MD", \
    ;                                      "GW RD (r="+sprintf("%4.2f",r_CTL_RD)+")",\
    ;                                      "FD RD  (r="+sprintf("%4.2f",r_FD_RD)+")",\
    ;                                      "GLEAM RD", \
    ;                                      "DOLCE RD"  /)
    textres@lgLabels                 = (/ \
                                        "FD MD  (r=0.28)",\
                                        "FD RD  (r=0.42)",\
                                        "GW MD (r=0.24)", \
                                        "GW RD (r=0.37)",\
                                        "GLEAM MD", \
                                        "GLEAM RD", \
                                        "DOLCE MD", \
                                        "DOLCE RD"  /)  
    textres@lgPerimOn                = False                   ; no perimeter
    textres@lgItemCount              = 8                       ; how many
    ; textres@lgItemOrder              = (/5, 7, 3, 1, 4, 6, 2, 0/)
    lineres                          = True
    lineres@lgLineLabelFontHeightF   = 0.015                   ; font height
    lineres@lgDashIndexes            = 0
    lineres@lgLineColors             = (/ "orange","gold","darkgreen","green3",  \
                                          "black","gray50","blue3","deepskyblue" \
                                       /)
    lineres@lgLineThicknesses        = 3.5                     ; line thicknesses

    plot = simple_legend(wks, plot, genres, lineres, textres)
    ; drawNDCGrid(wks)                                    ; This is for debugging purposes to help position the legend
    gsn_panel(wks,plot,(/1,1/),pres)
    ; frame(wks)
end if
end if

exit
end
