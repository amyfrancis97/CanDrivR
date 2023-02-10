#!/bin/bash

# Get conservation features
sbatch intersect.sh

# Get VEP features
sbatch getVEPCache.sh
sbatch queryVEPCache.sh

# Get dinucleotide properties
sbatch getDinucleotideProperties.sh

# Get amino acid properties
sbatch getAAProperties.sh

# Get DNA shape 
sbatch getDNAShape.sh
