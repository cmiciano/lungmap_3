#!/bin/sh
#PBS -q pdafm
#PBS -m abe
#PBS -M cmiciano@ucsd.edu
#PBS -V
#PBS -A epigen-group
#PBS -e /home/cmiciano/job_outs/Lung/lungmap_3/clust_multiome_filt.R.e
#PBS -o /home/cmiciano/job_outs/Lung/lungmap_3/clust_multiome_filt.R.o
#PBS -l walltime=05:00:00
#PBS -l nodes=1:ppn=4

# Navigate to oasis scratch

# activate cond env
source /home/cmiciano/miniconda3/etc/profile.d/conda.sh

conda activate ressen3

# clustering object with manually filtering of nCt_ATAC nFt_ATAC, TSS.enrichment and doublets removed 
# Input path

Script="/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/scripts/231109_clust_multiome.R"
Rscript $Script

