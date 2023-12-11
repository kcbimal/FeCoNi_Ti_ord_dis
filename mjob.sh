#!/bin/bash
#SBATCH -N 4
#SBATCH -p debug
#SBATCH --ntasks-per-node=40
#SBATCH --job-name=FeV_dis
#SBATCH -t 00:30:00
#SBATCH -o dis.out
#SBATCH -e dis.err
ulimit -s unlimited
module load intelmpi
export OMP_NUM_THREADS=1
echo "SLURM_NTASKS: " $SLURM_NTASKS

#python *.py
phonopy -p -s mesh.conf