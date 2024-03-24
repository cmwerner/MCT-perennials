#!/bin/bash -l

#SBATCH --account=coexistence
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cwerner5@uwyo.edu
#SBATCH --job-name=MCT-partitions-p

# Set the parameter combination to use and generate names of R scripts and log files
Rscript=mct-partitions-run-revision.R
LogFile=mct-partitions-run-revision.log

# Change to the relevant working directory
cd /project/coexistence/cwerner5

# Load R and MPI
module load gcc/12.2.0 openmpi/4.1.4 r-rmpi/0.7-1-ompi

R < $Rscript > $LogFile --no-save 
