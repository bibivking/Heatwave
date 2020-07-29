;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

begin

path1 = "/g/data/w35/mm3972/model/wrf/NUWRF/nuwrf_25km/LISWRF_configs/LIS_output"
path2 = "/g/data/w35/mm3972/model/cable/runs/my_version/run_Princeton/litter/outputs"

Soil_thickness  = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
Soil_thickness1 =(/0.005, 0.075, 0.154, 0.409, 1.085, 2.872/) ; LIS-CABLE
; region

A = -43.4377
B = -13.5853
C = 117.552
D = 153.598

;______________________________ Data _________________________________
year_s = 2003
year_e = 2008
year_sum = year_e-year_s+1
day_sum  = 6940

Var     = new((/2,20/),float)
Var     = 0.0

Var_name = (/\
             "Rain","Evap","TVeg","ESoil","ECanop","Qs","Qsb",\
             "Rnet","Qle","Qh","Qg","Tair","SoilMoist","SoilTemp",\
             "SoilMoist1","SoilMoist2","SoilMoist3","SoilMoist4", \
             "SoilMoist5","SoilMoist6"\
             /)

do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    ; read CABLE (HESS)
    filename1 = path2+"/cable_out_"+year+".nc"
    f1 = addfile (filename1,"r")

    ; read LIS-CABLE
    do month = 1, 12
      ;print("month = "+month)
      if month .lt. 10 then
        filename = path1+"/LIS.CABLE."+year+"0"+month+"0100.d01.nc"
      else
        filename = path1+"/LIS.CABLE."+year+month+"0100.d01.nc"
      end if

      ; print(filename)
      f = addfile (filename,"r")

      ; 0: Rain
      ;print(f->Rainf_f_inst)
      Var(0,0)     = Var(0,0) + avg(dim_sum_n_Wrap(f->Rainf_tavg,0))*3600*24
      ; 1: Evap
      Var(0,1)     = Var(0,1) + avg(dim_sum_n_Wrap(f->Evap_tavg,0))*3600*24
      ; 2: TVeg
      Var(0,2)     = Var(0,2) + avg(dim_sum_n_Wrap(f->TVeg_tavg,0))*3600*24
      ; 3: ESoil
      Var(0,3)     = Var(0,3) + avg(dim_sum_n_Wrap(f->ESoil_tavg,0))*3600*24
      ; 4: ECanop
      Var(0,4)     = Var(0,4) + avg(dim_sum_n_Wrap(f->ECanop_tavg,0))*3600*24
      ; 5: Qs
      Var(0,5)     = Var(0,5) + avg(dim_sum_n_Wrap(f->Qs_tavg,0))*3600*24
      ; 6: Qsb
      Var(0,6)     = Var(0,6) + avg(dim_sum_n_Wrap(f->Qsb_tavg,0))*3600*24
      ; 7: Rnet
      Var(0,7)     = Var(0,7) + avg(dim_avg_n_Wrap(f->Swnet_tavg,0))\
                              + avg(dim_avg_n_Wrap(f->Lwnet_tavg,0))
      ; 8: Qle
      Var(0,8)     = Var(0,8) + avg(dim_avg_n_Wrap(f->Qle_tavg,0))
      ; 9: Qh
      Var(0,9)     = Var(0,9) + avg(dim_avg_n_Wrap(f->Qh_tavg,0))
      ; 10: Qg
      Var(0,10)     = Var(0,10) + avg(dim_avg_n_Wrap(f->Qg_tavg,0))
      ; 11: Tair
      Var(0,11)     = Var(0,11) + avg(dim_avg_n_Wrap(f->Tair_f_inst,0))-273.15

      do soil = 0,5
        ; 12: SoilMoist
        Var(0,12)     = Var(0,12) + avg(dim_avg_n_Wrap(f->SoilMoist_inst(:,soil,:,:),0))\
                                   *Soil_thickness1(soil)/sum(Soil_thickness1)
        ; 13: SoilTemp
        Var(0,13)     = Var(0,13) + (avg(dim_avg_n_Wrap(f->SoilTemp_inst(:,soil,:,:),0))-273.15)\
                                   *Soil_thickness1(soil)/sum(Soil_thickness1)
        ; 14-19: SoilMoist1-6
        Var(0,14+soil)= Var(0,14+soil) + avg(dim_avg_n_Wrap(f->SoilMoist_inst(:,soil,:,:),0))
      end do ; soil
      delete(f)


      ; read CABLE (HESS)

      ; 0: Rain
      Var(1,0)     = Var(1,0) + avg(f1->Rainf(month-1,{A:B},{C:D})*dom(month-1))*3600*24
      ; 1: Evap
      Var(1,1)     = Var(1,1) + avg(f1->Evap(month-1,{A:B},{C:D})*dom(month-1))*3600*24
      ; 2: TVeg
      Var(1,2)     = Var(1,2) + avg(f1->TVeg(month-1,{A:B},{C:D})*dom(month-1))*3600*24
      ; 3: ESoil
      Var(1,3)     = Var(1,3) + avg(f1->ESoil(month-1,{A:B},{C:D})*dom(month-1))*3600*24
      ; 4: ECanop
      Var(1,4)     = Var(1,4) + avg(f1->ECanop(month-1,{A:B},{C:D})*dom(month-1))*3600*24
      ; 5: Qs
      Var(1,5)     = Var(1,5) + avg(f1->Qs(month-1,{A:B},{C:D})*dom(month-1))*3600*24
      ; 6: Qsb
      Var(1,6)     = Var(1,6) + avg(f1->Qsb(month-1,{A:B},{C:D})*dom(month-1))*3600*24


    end do ; month

    ; 7: Rnet
    Var(1,7)     = Var(1,7) + avg(dim_avg_n_Wrap(f1->Rnet(:,{A:B},{C:D}),0))
    ; 8: Qle
    Var(1,8)     = Var(1,8) + avg(dim_avg_n_Wrap(f1->Qle(:,{A:B},{C:D}),0))
    ; 9: Qh
    Var(1,9)     = Var(1,9) + avg(dim_avg_n_Wrap(f1->Qh(:,{A:B},{C:D}),0))
    ; 10: Qg
    Var(1,10)     = Var(1,10) + avg(dim_avg_n_Wrap(f1->Qg(:,{A:B},{C:D}),0))
    ; 11: Tair
    Var(1,11)     = Var(1,11) + avg(dim_avg_n_Wrap(f1->Tair(:,{A:B},{C:D}),0))-273.15

    do soil = 0,5
      ; 12: SoilMoist
      Var(1,12)     = Var(1,12) + avg(dim_avg_n_Wrap(f1->SoilMoist(:,soil,{A:B},{C:D}),0))\
                                 *Soil_thickness(soil)/sum(Soil_thickness)
      ; 13: SoilTemp
      Var(1,13)     = Var(1,13) + (avg(dim_avg_n_Wrap(f1->SoilTemp(:,soil,{A:B},{C:D}),0))-273.15)\
                                 *Soil_thickness(soil)/sum(Soil_thickness)
      ; 14-19: SoilMoist1-6
      Var(1,14+soil)= Var(1,14+soil) + avg(dim_avg_n_Wrap(f1->SoilMoist(:,soil,{A:B},{C:D}),0))
    end do ; soil
    delete(f1)
end do ; year
print(Var)
Var(0,0:6) = Var(0,0:6)/year_sum
Var(0,7:) = Var(0,7:)/year_sum/12
Var(1,0:6) = Var(1,0:6)/year_sum
Var(1,7:) = Var(1,7:)/year_sum

do var_num = 0, dimsizes(Var_name)-1
  print(Var_name(var_num)+ " "+ Var(0,var_num)+" "+ Var(1,var_num))
end do

end
