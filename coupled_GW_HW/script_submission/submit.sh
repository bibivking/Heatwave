#!/bin/bash

#PBS -m ae
#PBS -P w35
#PBS -q express
#PBS -l walltime=4:00:00
#PBS -l mem=190GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w35

source activate science
cd /g/data/w35/mm3972/scripts/Heatwave/coupled_GW_HW
python Fig1_spatial_map_wrf_domain_pft.py
