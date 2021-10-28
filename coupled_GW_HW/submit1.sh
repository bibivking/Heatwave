#!/bin/bash

#PBS -m ae
#PBS -P w35
#PBS -q normalbw
#PBS -l walltime=0:10:00
#PBS -l mem=256GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w35

source activate science
cd /g/data/w35/mm3972/scripts/Heatwave/coupled_GW_HW
python spatial_map_weather_analysis.py
