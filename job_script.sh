#!/bin/bash
#SBATCH --job-name=cyounes-job  # Job name
#SBATCH --account=randommatrix  # Account to charge
#SBATCH --partition=compute     # Partition to use
#SBATCH --nodes=4               # Number of nodes
#SBATCH --ntasks=4              # One task per node
#SBATCH --cpus-per-task=40      # 40 CPUs per node
#SBATCH --mem=175G              # Memory allocation
#SBATCH --time=10:00:00         # Maximum runtime
#SBATCH --output=cyounes-job-%j.out  # Standard output log
#SBATCH --error=cyounes-job-%j.err   # Standard error log

# Load necessary modules (adjust as needed)
module load elbert/julia/1.10.2/1.10.2

# Create a list of allocated nodes
scontrol show hostnames > hosts.txt
cat hosts.txt  # Print node list for debugging

# Set Julia environment
export JULIA_NUM_THREADS=40  # Each process will use 40 threads

# Run Julia with Distributed Computing
nohup srun julia -p 160 Example1.jl > output.log 2>&1 &