#!/bin/bash
#PBS -q condo
#PBS -m abe
#PBS -M micianoch@gmail.com
#PBS -V
#PBS -A epigen-group
#PBS -e /home/cmiciano/job_outs/Lung/lungmap_3/df_atac.R.e
#PBS -o /home/cmiciano/job_outs/Lung/lungmap_3/df_atac.R.o
#PBS -l walltime=00:40:00
#PBS -l nodes=1:ppn=4

# generally 45 min walltime

# include these 2 lines since bash commands aren't recognized outside my base environment (subshells)
source /home/cmiciano/miniconda3/etc/profile.d/conda.sh
conda activate ressen3

# input: txt file space delimited
# <location of filtered .h5" <min_nFeature_RNA> <max_nFeature_RNA> <percent.mt> <condition> <age at donation> <age> <age units> <sex> <race>
# JB_224/JB_224_filtered.h5 100 4000 8 BPD_Chronic Child 15 Months Male White

INPUT_FILE="/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/scripts/240104_LungMap3_PairedTag_Multiome_ATAC_metadata_sp_fi.tsv"
SAMPLE_ROW=`cat $INPUT_FILE | sed -n ${PBS_ARRAYID}p`


# Using the absolute paths to .h5s
SAMP_FP=$(echo $SAMPLE_ROW | cut -f1 -d' ')

# Outdir for objects and plots
OUT_PATH="/oasis/tscc/scratch/cmiciano/Lung/lungmap_3/02_DoubletFinder/"

[[ -d $OUT_PATH ]] || mkdir $OUT_PATH


# sample name (library ID ie JB_XX)
LIB_RNA_ID=$(echo $SAMPLE_ROW | cut -f2 -d' ')

# Sample ID according to collaborator (adam's id)
SAMP_ID=$(echo $SAMPLE_ROW | cut -f3 -d' ')

MODALITY="multiome"

LIB_ATAC_ID=$(echo $SAMPLE_ROW | cut -f4 -d' ')

LIB_RNA_ATAC_ID=$(echo $SAMPLE_ROW | cut -f5 -d' ')

FRAG_PATH=$(echo $SAMPLE_ROW | cut -f6 -d' ')

echo "samp_fp: " $SAMP_FP
echo "lib rna id: " $LIB_RNA_ID 
echo "samp_id: " $SAMP_ID
echo "modality" $MODALITY

echo "out path: " $OUT_PATH

echo "lib atac id: " $LIB_ATAC_ID 
echo "lib rna atac id: " $LIB_RNA_ATAC_ID 

echo "frag_path" $FRAG_PATH
#echo "ndims" $NDIMS

# #Rscript scRNA_prep.R H5 $H5 SP "mm10" LIB "JB_922_1_2" #OUT $OUT
# Rscript scRNA_prep.R -h5 $H5 -l "JB_922_1_2" -s "mm10" -o $OUT

Script='/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/scripts/231026_doubletFinder.R'

Rscript $Script --inp_h5 $SAMP_FP --lib_rna_id $LIB_RNA_ID --samp_id $SAMP_ID --modality $MODALITY --out $OUT_PATH --lib_atac_id $LIB_ATAC_ID --lib_rna_atac_id $LIB_RNA_ATAC_ID --fragpath $FRAG_PATH
#                  # 1         #2      #3           #4         #5       #6         #7             #8
#Rscript $Script $SAMP_FP $OUT_PATH $LIB_RNA_ID $SAMP_ID $TISSUE $CAGE $SEX $AGE $TREATMENT $TREATMENT_A_B $AGE_DAYS $BIRTHDATE $LIB_ATAC_ID $LIB_RNA_ATAC_ID $MODALITY $FRAG_PATH $NDIMS
                    #1        #2        #3      #4        #5      #6    #7   #8        #9              #10    #11      #12         #13        #14              #15      #16        #17
# Copy over the output files for the single sample to ps-epigen
for DIR in $OUT_PATH*;do
        echo $DIR;
        BN=$(basename $DIR);
        echo $BN;
        FOUT="/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/02_DoubletFinder/"
        [[ -d $FOUT ]] || mkdir $FOUT
        FDOUT=$FOUT$BN
        [[ -d $FDOUT ]] || mkdir $FDOUT
        echo $FDOUT
        cp -R $DIR/$LIB_RNA_ATAC_ID* $FDOUT

done
