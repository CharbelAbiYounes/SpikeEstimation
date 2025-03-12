#!/bin/bash
#SBATCH --job-name=cyounes2-job  # Job name
#SBATCH --account=amath  # Account to charge
#SBATCH --partition=cpu-g2     # Partition to use
#SBATCH --nodes=1               # Request 1 node
#SBATCH --cpus-per-task=80      # Request CPUs
#SBATCH --mem=1400G              # Memory allocation
#SBATCH --time=24:00:00         # Maximum runtime
#SBATCH --output=cyounes-Distjob-%j.out  # Standard output log
#SBATCH --error=cyounes-Distjob-%j.err   # Standard error log

module load elbert/julia/1.10.2/1.10.2

# scontrol show hostnames > hosts.txt
# cat hosts.txt  

export JULIA_NUM_THREADS=80

srun taskset -c 0-79 nice -n 0 julia Example1.jl > output.log 2>&1