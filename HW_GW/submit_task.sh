#!/bin/bash

#PBS -m ae
#PBS -M mengyuan.mu@unsw.edu.au
#PBS -P w35
#PBS -q normal
#PBS -l walltime=05:00:00
#PBS -l mem=40GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w35

module load ncl/6.6.2

cd /g/data/w35/mm3972/scripts/Heatwave/HW_GW

ncl boxplot_HESS_HW_2000-2019.ncl
