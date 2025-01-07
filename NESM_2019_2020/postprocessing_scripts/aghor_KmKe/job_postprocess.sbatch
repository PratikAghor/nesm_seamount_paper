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
WORKDIR=/storage/home/hcoda1/9/paghor3/scratch/work/NESM_2019_2020/postprocessing_scripts/aghor_KmKe
#WORKDIR=./
cd $WORKDIR

# save annual mean
# matlab -nodesktop -nosplash -nodisplay < save_annual_mean.m

# save hrs
# matlab -nodesktop -nosplash -nodisplay < get_hrs.m

# save vrs
# matlab -nodesktop -nosplash -nodisplay < get_vrs.m

# save hrs 3d
# matlab -nodesktop -nosplash -nodisplay < save_hrs_3d.m

# save vrs 3d
# matlab -nodesktop -nosplash -nodisplay < save_vrs_3d.m

# take vertically integrated 2d hslice of hrs and vrs, requires hrs_3d, vrs_3d files
# matlab -nodesktop -nosplash -nodisplay < save_hrs_2d.m
# matlab -nodesktop -nosplash -nodisplay < save_vrs_2d.m

# get 3d KmKe timeseries
matlab -nodesktop -nosplash -nodisplay < save_KmKe_3d.m

# get 2d KmKe timeseries from the 3d timeseries
matlab -nodesktop -nosplash -nodisplay < save_KmKe_2d.m                                                                                                                                           
