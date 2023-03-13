#!/bin/bash

#sbatch 1_getConservationFeatures.sh coding
#sbatch 1_getConservationFeatures.sh non-coding

#sbatch 2_getVEPCache.sh

#sbatch 3_queryVEPCache.sh coding
#sbatch 3_queryVEPCache.sh non-coding

sbatch 4_getDinucleotideProperties.sh coding
sbatch 4_getDinucleotideProperties.sh non-coding

#sbatch 5_getDNAShape.sh coding
#sbatch 5_getDNAShape.sh non-coding

#sbatch 6_getGC_GpG.sh coding
#sbatch 6_getGC_GpG.sh non-coding

#sbatch 7_getkernel.sh coding
#sbatch 7_getkernel.sh non-coding
