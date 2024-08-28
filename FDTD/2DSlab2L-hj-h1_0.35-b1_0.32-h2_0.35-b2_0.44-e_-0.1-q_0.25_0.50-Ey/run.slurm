#!/bin/bash
#SBATCH --qos=long
#SBATCH --job-name=FDTD-e_-0.1-Ey
#SBATCH --output=FDTD-e_-0.1.o
#SBATCH --error=FDTD-e_-0.1.e
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --mail-user=ngocduc.le@dipc.org
#SBATCH --time=2-00:00:00 

module load Python
conda activate /scratch/lengocduc/.conda/envs/pmp

set OMP_NUM_THREAS= $SLURM_NTASKS

time mpirun -np $SLURM_NTASKS python 2DSlab2L-hj-qd-6.py >& output

conda deactivate 
