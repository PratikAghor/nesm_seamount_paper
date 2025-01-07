#!/bin/bash
#SBATCH -Jpostprocess
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=24G 
#SBATCH -t24:00:00
#SBATCH -oReport-%j.out                  # Combined output and error messages file
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
# WORKDIR=/storage/home/hcoda1/9/paghor3/scratch/work/SEAMOUNT_2019_2020/postprocessing_scripts/
WORKDIR=./
cd $WORKDIR
# save pv_3d_avg
# matlab -nodesktop -nosplash -nodisplay < save_pv_3d_avg.m

# save pv_3d_t -> timeseries of 3d epv
matlab -nodesktop -nosplash -nodisplay < save_pv_3d_t.m


