#!/bin/bash
#SBATCH -J uniform                   # A single job name for the array
#SBATCH -n 12                      # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH -p serial_requeue   # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 0-1:10                  # Maximum execution time (D-HH:MM)
#SBATCH -o out/1router_t5_uniform.out           # Standard output
#SBATCH -e dump/1router_t5_uniform.err          # Standard error

module load R/3.2.0-fasrc01
Rscript wonlee_1router_t5_uniform.R
