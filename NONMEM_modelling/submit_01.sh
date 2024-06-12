#!/bin/bash -l

#SBATCH --cluster=wice
#SBATCH --account=lp_dreesenlab
#SBATCH --nodes=1
#SBATCH --partition=batch_sapphirerapids
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=96
#SBATCH --time=02:00:00

module load NONMEM/7.5.0-GCC-10.3.0-MPICH-3.4.2
module load PsN/5.3.0-foss-2021a

cd $SLURM_SUBMIT_DIR
execute -nodes=1 -threads=1 -nm_version=nm750 -parafile=mpilinux8.pnm run_onco_007.mod