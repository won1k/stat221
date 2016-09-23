#!/bin/bash
#SBATCH -J task3                   # A single job name for the array
#SBATCH -n 12                      # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH -p serial_requeue          # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 0-1:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o out/task3.out           # Standard output
#SBATCH -e dump/task3.err          # Standard error

