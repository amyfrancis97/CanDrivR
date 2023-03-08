#!/bin/bash
#SBATCH --job-name=getDinucleotideProperties
#SBATCH --partition=cnu
#SBATCH --mem=100G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

# Download dinucleotide property table from https://diprodb.fli-leibniz.de/ShowTable.php
module load lang/r

Rscript dinucleotideProperties.R
