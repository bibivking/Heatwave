#!/bin/bash

#PBS -m ae
#PBS -M mengyuan.mu@unsw.edu.au
#PBS -P w35
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l mem=180GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w35

module load ncl/6.6.2

cd /g/data/w35/mm3972/scripts/Heatwave/HW_GW

#ncl scatter_HESS_HW_deltaT_WTD.ncl

ncl time_accumulate_HESS_P-E.ncl
