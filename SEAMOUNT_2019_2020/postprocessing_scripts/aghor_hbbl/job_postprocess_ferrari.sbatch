#!/bin/bash
#SBATCH -Jhbbl_ferrari
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=96G 
#SBATCH -t96:00:00
#SBATCH -oReport-%j.out                  # Combined output and error messages file
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
# WORKDIR=/storage/home/hcoda1/9/paghor3/scratch/work/SEAMOUNT_2019_2020/postprocessing_scripts/aghor_KmKe
WORKDIR=./
cd $WORKDIR

# save annual mean
# matlab -nodesktop -nosplash -nodisplay < save_rho_w_annual_mean.m

# save hbbl_t
# matlab -nodesktop -nosplash -nodisplay < save_hbbl_t.m

# save Ec_2d_avg
# matlab -nodesktop -nosplash -nodisplay < save_Ec_2d.m

# save kv_t_vslices
# matlab -nodesktop -nosplash -nodisplay < save_kv_t_vslice_ferrari.m

matlab -nodesktop -nosplash -nodisplay < save_kv_t_vslice_laurent_ferrari.m
