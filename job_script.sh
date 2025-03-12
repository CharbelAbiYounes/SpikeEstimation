#!/bin/bash
#SBATCH --job-name=cyounes-job  # Job name
#SBATCH --account=randommatrix  # Account to charge
#SBATCH --partition=compute     # Partition to use
#SBATCH --nodes=4               # Number of nodes
#SBATCH --ntasks=4              # One task per node
#SBATCH --cpus-per-task=8      # 40 CPUs per node
#SBATCH --mem=175G              # Memory allocation
#SBATCH --time=48:00:00         # Maximum runtime
#SBATCH --output=cyounes-job-%j.out  # Standard output log
#SBATCH --error=cyounes-job-%j.err   # Standard error log

# Load the Julia module
module load elbert/julia/1.10.2/1.10.2

scontrol show hostnames > hosts.txt
cat hosts.txt  

# Set the number of Julia threads and assign CPUs to the task
export JULIA_NUM_THREADS=40

# Run the Julia script with taskset to bind the process to specific CPUs
srun taskset -c 0-39 nice -n 0 julia Example1.jl > output.log 2>&1