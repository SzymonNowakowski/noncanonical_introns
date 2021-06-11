#!/bin/bash -l
#SBATCH -J gpu_test
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1 
#SBATCH -A GS82-12
#SBATCH --time=02:00:00 


ln -s ~/introns /scratch/$SLURM_JOB_ID/introns 


singularity exec --nv  --bind /scratch/$SLURM_JOB_ID:/scratch docker://akalinow/tensorflow-gpu:latest /scratch/introns/noncanonical_introns/runme.sh >& ~/output-$SLURM_JOB_ID

