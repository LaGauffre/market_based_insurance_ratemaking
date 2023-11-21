#!/bin/bash


# number of cores
#SBATCH -n 128
# Walltime
#SBATCH -t 0:10:00
# send a mail at the end of the exec
#SBATCH --mail-type=END
#SBATCH --mail-user=goffard@unistra.fr

R CMD BATCH poisson_lognormal.R

