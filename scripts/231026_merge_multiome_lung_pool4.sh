#!/bin/bash
#PBS -q home-epigen
#PBS -m abe
#PBS -M cmiciano@ucsd.edu
#PBS -V
#PBS -A epigen-group
#PBS -e /home/cmiciano/job_outs/Lung/lungmap_3/multiome_merg_lungpool4.R.e
#PBS -o /home/cmiciano/job_outs/Lung/lungmap_3/multiome_merg_lungpool4.R.o
#PBS -l walltime=05:00:00
#PBS -l nodes=1:ppn=28


# liver seems to fail with 16ppn in hotel
# merging multiome samples 
# atac takes 15 min per sample to quantify
source /home/cmiciano/miniconda3/etc/profile.d/conda.sh
conda activate ressen3


# testing on list of lung objects
# lung 30 pool 4 samples
RDS_PATH="/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/01_DoubletFinder/df_objs/"

OUT="/oasis/tscc/scratch/cmiciano/Lung/lungmap_3/03_multiome_merge/" ##don't forget the slash

PROJ="231026_merged_multiome_lung4"

Script="/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/scripts/231026_merge_multiome.R"


echo "filepaths: " $RDS_PATH
echo "outdir: "  $OUT
echo "project: " $PROJ

[[ -d $OUT ]] || mkdir $OUT

Rscript $Script $RDS_PATH $OUT $PROJ


## Copy outputs from scratch to ps-epigen
FOUT="/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/"
cp -R $OUT $FOUT
