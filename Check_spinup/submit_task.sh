#!/bin/bash

#PBS -m ae
#PBS -M mengyuan.mu@unsw.edu.au
#PBS -P w35
#PBS -q normal
#PBS -l walltime=03:00:00
#PBS -l mem=20GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w35

module load ncl/6.6.2

cd /g/data/w35/mm3972/scripts/Heatwave/Check_spinup
ncl spinup_month-daily_CABLE-HESS.ncl
