#!/bin/bash
#SBATCH -Jreynoldsstresses
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=24G 
#SBATCH -t96:00:00
#SBATCH -oReport-%j.out                  # Combined output and error messages file
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
WORKDIR=/storage/home/hcoda1/9/paghor3/scratch/work/SEAMOUNT_2019_2020/postprocessing_scripts/aghor_KmKe
#WORKDIR=./
cd $WORKDIR

# save hrs
matlab -nodesktop -nosplash -nodisplay < get_hrs.m

# save vrs
matlab -nodesktop -nosplash -nodisplay < get_vrs.m
