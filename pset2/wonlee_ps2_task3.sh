#!/bin/bash
#SBATCH -J task3                   # A single job name for the array
#SBATCH -n 12                      # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH -p stats   # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 0-1:10                  # Maximum execution time (D-HH:MM)
#SBATCH -o out/task3_%A_%a.out           # Standard output
#SBATCH -e dump/task3_%A_%a.err          # Standard error

module load R/3.2.0-fasrc01
Rscript wonlee_ps2_task3.R "${SLURM_ARRAY_TASK_ID}"
