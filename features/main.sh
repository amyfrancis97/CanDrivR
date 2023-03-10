#!/bin/bash

#sbatch 1_getConservationFeatures.sh coding
#sbatch 1_getConservationFeatures.sh non-coding

#sbatch 2_getVEPCache.sh

#sbatch 3_queryVEPCache.sh coding
#sbatch 3_queryVEPCache.sh non-coding

#sbatch 4_getDinucleotideProperties.sh coding
sbatch 4_getDinucleotideProperties.sh non-coding
