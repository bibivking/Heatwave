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

cd /g/data/w35/mm3972/scripts/Heatwave/HW_GW/run_multi_scripts

ncl pdf_calc_HESS_HW_Qh.ncl
