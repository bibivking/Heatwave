#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q express
#PBS -l walltime=4:00:00
#PBS -l mem=190GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97+gdata/hh5

module load conda/analysis3-22.04
cd /g/data/w97/mm3972/scripts/Heatwave/coupled_GW_HW
python Fig4_profile_wrf_var_transect.py
