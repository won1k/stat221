#!/bin/bash
#SBATCH -J informative                   # A single job name for the array
#SBATCH -n 12                      # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH -p serial_requeue   # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 0-1:10                  # Maximum execution time (D-HH:MM)
#SBATCH -o out/1router_all_uniform_%A_%a.out           # Standard output
#SBATCH -e dump/1router_all_uniform_%A_%a.err          # Standard error

module load R/3.2.0-fasrc01
Rscript wonlee_1router_all_times_uniform.R $SLURM_ARRAY_TASK_ID
