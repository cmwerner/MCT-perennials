#!/bin/bash -l

#SBATCH --account=coexistence
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cwerner5@uwyo.edu
#SBATCH --job-name=MCT-partitions-p

# Set the parameter combination to use and generate names of R scripts and log files
Rscript=mct-partitions-run-p.R
LogFile=mct-partitions-run-p.log

# Change to the relevant working directory
cd /project/coexistence/cwerner5

# Load R and MPI
module load gcc/7.3.0 r/3.5.3 swset/2018.05 openmpi/3.1.0 r-rmpi/0.6-9-r353-py27

R < $Rscript > $LogFile --no-save 
