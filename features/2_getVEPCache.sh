#!/bin/bash
#SBATCH --job-name=getVEPCache
#SBATCH --partition=gpu
#SBATCH --mem=50G
#SBATCH --time=1-00:00:0
#SBATCH --chdir=/bp1/mrcieu1/users/uw20204/paper1/features
#SBATCH --account=sscm013903

# Downloads the cache for VEP
curl -L https://cpanmin.us | perl - App::cpanminus
cpanm Archive::Zip
cpanm DBD::mysql
module load apps/bayestraits apps/bcftools apps/samtools apps/tabix lib/htslib

# There are a large number of package requirements that must be met before the below script can be carried out
# Please follow package dependencies at this link:https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer
git clone https://github.com/Ensembl/ensembl-vep
cd ensembl-vep
perl INSTALL.pl # Download the GRCh38 genome when prompted

