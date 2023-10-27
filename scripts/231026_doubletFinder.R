suppressMessages(library(rmarkdown))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(rhdf5))
suppressMessages(library(ggplot2))
suppressMessages(library(sctransform))
suppressMessages(library(DoubletFinder))
suppressMessages(library(DT))
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
suppressMessages(library(hdf5r))
library(Signac)
library(argparse)


# What does this do?
# Takes the .h5 file from cellbender and run's doubletFinder on it

# Why are we doing it?
# Mark doublets
# Adds condition and sample metadata

# How does it select doublets?
# it takes the number of cells present after filtering, 
# 

# How does it filter low quality cells?
# it calculates the quantiles for the percent mito and the nFeatures
# #shave bottom 5% and top 5% of nFeatures and top 5% of mito
# generally 80-85% cells remain in each sample after filtering

# How to run? qsub doubletFinder -t <index of sample 1> <index of last sample>

# How long does it take?
# Around between 15-30 min


## Creating seurat object 

# 
find_nexp <- function(inp_cells) {

 rec_cells <- round_any(inp_cells, 1000, f=floor)

  # If number of cells inputted exceed table, manually set percent
  if(inp_cells > 10000) {
     chos_pt <- 0.076
     rec_cells <- 10000
  } else if(inp_cells < 1000){
     chos_pt <- 0.004
     rec_cells <- 500
  }
  else {
     mtable <- read.delim("/projects/ps-epigen/users/cmiciano/SenNet_Multiome/scripts/multiplet_rate.txt", sep=" ", header = T, stringsAsFactors = F)
     chos_pt <- mtable[which(mtable$recov_cells == rec_cells),c("percent")]
  }

     #Calculate percentage based on 10X multiplet rate table
     print(paste0("Chosen percent: ", chos_pt))
     calc_pt <- (inp_cells*chos_pt)/rec_cells
     print(paste0("Calculated percentage: ", calc_pt))

     #nExp for doublet finder param
     print(paste0("Number of cells after filtering: ", inp_cells))
     calc_nexp <- round(calc_pt*inp_cells) #multiply by number of cells after filtering instead of the estimated number of cells inputted as cellbender params
     print(paste0("Calculated nExp: ", calc_nexp))
     return(calc_nexp)
}




# $SAMP_FP $NFEAT_RNA_MIN $NFEAT_RNA_MAX $PERCENT_MT $NUM_DIMS $PN $OUT_PATH
# args
#args = commandArgs(trailingOnly=TRUE)
#inp_h5 <- args[1]
#out <- args[2]
#lib_rna_id <- args[3]
#samp_id <- args[4]
#tissue <- args[5]
#cage <- args[6]
#sex <- args[7]
#age <- args[8]
#treatment <- args[9]
#treatment_a_b <- args[10]
#age_days <- args[11]
#birthdate <- args[12]
#lib_atac_id <- args[13]
#lib_rna_atac_id <- args[14]
#modality <- args[15]
#fragpath <- args[16]
#ndims <- args[17]
# library_id


#print(paste("h5 ", inp_h5))
#print(paste("outpath ", out))
#print(paste("lib_rna_id ", lib_rna_id))
#print(paste("samp_id ", samp_id))
#print(paste("tissue ", tissue))
#print(paste("cage ", cage))
#print(paste("sex ", sex))
#print(paste("age ", age))
#print(paste("treatment ", treatment))
#print(paste("treatment_a_b ", treatment_a_b))
#print(paste("age_days ", age_days))
#print(paste("birthdate ", birthdate))
#print(paste("lib_atac_id ", lib_atac_id))
#print(paste("lib_rna_atac_id ", lib_rna_atac_id))
#print(paste("modality ", modality))
#print(paste("fragpath ", fragpath))
#print(paste("ndims ", ndims))


#Rscript $Script $SAMP_FP $LIB_RNA_ID $SAMP_ID $MODALITY $OUT_PATH $LIB_ATAC_ID $LIB_RNA_ATAC_ID $FRAGPATH
#                  # 1         #2      #3           #4         #5       #6         #7             #8
# create parser object
parser <- ArgumentParser(description = "Create seurat object from 10x h5 and collect QC metrics")

required_arg_group = parser$add_argument_group('flagged required arguments',
        'the script will fail if these args are not included')

# These args belong to my group of required flagged arguments
required_arg_group$add_argument("-h5", "--inp_h5", required = TRUE,
    help="path to dataset h5 file for object creation")

#required_arg_group$add_argument("-s","--species", required = TRUE,
#    help="species of origin for the dataset (mm10 or hg38)")

required_arg_group$add_argument("-lr","--lib_rna_id", required = TRUE,
     help="library rna_id of the sample (JB_232)")

required_arg_group$add_argument("-s","--samp_id", required = TRUE,
     help="library rna_id of the sample (JB_232)")

required_arg_group$add_argument("-m","--modality", required = TRUE,
     help="modality of RNA or multiome")

required_arg_group$add_argument("-o","--out", required = TRUE,
    help="output path for plots and objects")

parser$add_argument("-skipf","--skip_filter", action = 'store_true',
    help="determines whether to skip filter")

parser$add_argument("-skip","--skip_df", action = 'store_true',
    help="determines whether to skip doubletfinder and only filter")

parser$add_argument("-la","--lib_atac_id", required = FALSE,
     help="library atac_id of the sample (JB_231)")

parser$add_argument("-lra","--lib_rna_atac_id", required = FALSE,
     help="combined library rna_id and atac id of the sample (JB_232_1_2_JB_231_1_2)")

parser$add_argument("-f","--fragpath", required = FALSE,
     help="fragment path to atac for multiome objs")
args <- parser$parse_args()


print("args")
print(args$inp_h5)
print(args$lib_rna_id)
print(args$samp_id)
print(args$modality)
print(args$out)
print(args$lib_atac_id)
print(args$lib_rna_atac_id)
print(args$fragpath)
 

# Log inputs
CB_file <- args$inp_h5
CB_data <- Read10X_h5(filename = CB_file, use.names = TRUE)

if (args$modality == "multiome") {
  # get gene annotations for hg38
  #annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, verbose = F)
  #seqlevelsStyle(annotation) <- "UCSC"

  sobj_raw <- CreateSeuratObject(
   counts = CB_data$`Gene Expression`,
   assay = "RNA", project = args$lib_rna_atac_id
  )


  # create ATAC assay and add it to the object
  sobj_raw[["ATAC"]] <- CreateChromatinAssay(
    counts = CB_data$Peaks,
    sep = c(":", "-"),
    fragments = args$fragpath,
    min.features = -1, #important or else number of cells will be different
    #annotation = annotation
   )

   # Calculating frip
   total_fragments <- CountFragments(args$fragpath, cells = colnames(sobj_raw))

   frg <- total_fragments[ , c("CB","frequency_count")]
   rownames(frg) <- frg$CB
   frg$CB <- NULL
   sobj_raw <- AddMetaData(
   	 object = sobj_raw,
  	 metadata = frg,
         col.name = 'fragments_freq_count'
   )
   sobj_raw <- FRiP(
         object = sobj_raw,
         assay = 'ATAC',
         total.fragments = 'fragments_freq_count'
   )	
   


} else {

  sobj_raw <- CreateSeuratObject(counts = CB_data, project = args$lib_rna_id)

}

# set default assay to rna for proper percent.mt converation
DefaultAssay(sobj_raw) <- 'RNA'
## Filtering nFeatures and percent.mt

sobj_raw[["percent.mt"]] <- PercentageFeatureSet(sobj_raw, pattern = "^MT-")

# adding ATAC features



print("Before Filtering:")
print(paste("Number Nuclei Before Filtering: ", nrow(sobj_raw@meta.data), sep = " "))
bf_filt <- VlnPlot(sobj_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Filtering by quantile

if (args$skip_filter == FALSE) { 

quants <- quantile(sobj_raw@meta.data$nFeature_RNA, c(0.05, 0.95))
quants

pct_five <- quants[1]
nFeatRNA_min <- ceiling(pct_five)
print(paste("Keeping above nFeat_RNA: ", nFeatRNA_min, sep = " "))

pct_ninefive <- quants[2]
nFeatRNA_max <- ceiling(pct_ninefive)
print(paste("Keeping below nFeat_RNA: ", nFeatRNA_max, sep = " "))

sub_feat <- subset(sobj_raw, nFeature_RNA > nFeatRNA_min & nFeature_RNA < nFeatRNA_max) #shave bottom 5% and top 5%

quants_mito <- quantile(sub_feat@meta.data$percent.mt, c(0.95)) #shave top 5%
pct_mito_ninefive <- quants_mito[1]
percent_mt <- ceiling(pct_mito_ninefive)
if (percent_mt == 0) {
    percent_mt <- 1
}

#print(paste("Keeping below percent_mt: ", percent_mt, sep = " "))
#filt_mito_feat <- subset(sub_feat, percent.mt < percent_mt) 

# for now not filtering at all
sobj <- sobj_raw

#sobj <- subset(sobj_raw, subset = nFeature_RNA > nFeatRNA_min & nFeature_RNA < nFeatRNA_max & percent.mt < percent_mt)
print("After Filtering:")
print(paste("Number Nuclei After Filtering: ", nrow(sobj@meta.data), sep = " "))

raw_cell_cts <-  nrow(sobj_raw@meta.data)
filt_cell_cts <- nrow(sobj@meta.data)

pct_cells_kept <- round(filt_cell_cts / raw_cell_cts, 4) * 100
af_filt <- VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
af_filt <- af_filt + plot_annotation(title=paste0("ncells after filt: ", "n=", filt_cell_cts, " orig=", raw_cell_cts,"\n", "pct=", pct_cells_kept),
             theme = theme(plot.title = element_text(hjust = 0.5)))

dir.create(path= paste0(args$out), showWarnings= FALSE)

## Create filtered plots
print("Creating filtered plot")
dir.create(path= paste0(args$out,"filt_vlnplots"), showWarnings= FALSE)
pdf(file=paste0(args$out,"filt_vlnplots/",args$lib_rna_atac_id, "_vlnplot_filt.pdf"), title= paste0(args$lib_rna_atac_id,"_vlnplot_filt"), width= 12, height= 10, compress=FALSE)
bf_filt
af_filt
dev.off()

}

# Processing
# jb185
print("Processing...")
sobj <- SCTransform(sobj, verbose = FALSE)
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj), npcs = 50, verbose = FALSE)
sobj <- RunUMAP(sobj, dims = 1:50)


# Finding pK

# Don't run doubletFinders

if (args$skip_df == FALSE) {

print("Running DoubletFinder")
#print("Finding pK...")
sweep.sobj <- paramSweep_v3(sobj, PCs = 1:15, sct = TRUE)
sweep.stats.sobj <- summarizeSweep(sweep.sobj, GT = FALSE)
bcmvn.sobj <- find.pK(sweep.stats.sobj)
datatable(bcmvn.sobj, rownames = F)

# choosing maxima of BCmetric
pK <- as.numeric(as.character(bcmvn.sobj$pK))
top_pK <- pK[which(bcmvn.sobj$BCmetric %in% max(bcmvn.sobj$BCmetric))]
print(paste("Potential pK:", top_pK, sep = " "))


## You will see this "bcmvn.sobj" is a dataframe with 5 columns, which are: ParamID, pK, MeanBC,  VarBC, BCmetric

## Note that the "pK" column could be a factor, so converted to character and then numeric might be necessary.
dir.create(path= paste0(args$out,"pK_plots"), showWarnings= FALSE)
pdf(file=paste0(args$out,"pK_plots/",args$lib_rna_atac_id, "_pKplot.pdf"), title= paste0(args$lib_rna_atac_id,"_pKplot"), width= 12, height= 10, compress=FALSE)
plot(x = as.numeric(as.character(bcmvn.sobj$pK)), y = as.numeric(bcmvn.sobj$BCmetric), pch = 16, col = "#41b6c4", cex = 0.75, xlab="pK", ylab="BCmetric")
lines(x = as.numeric(as.character(bcmvn.sobj$pK)), y = bcmvn.sobj$BCmetric, col = "#41b6c4")
title(main = paste0(args$lib_rna_atac_id," ", "Top pK: ", top_pK)) 
abline(v=top_pK,lwd=2,col='red',lty=2)
dev.off()


# Calculating nEXP
# jb185 (est # cells - 8644)

print("Calculating nExp...")
cells_aft_filt <- nrow(sobj@meta.data)
calc_nExp <- find_nexp(cells_aft_filt)
nExp_poi_sobj <- round(calc_nExp*nrow(sobj@meta.data))
print(paste0("Calculated nExp: ", calc_nExp))

sobj <- doubletFinder_v3(sobj, PCs = 1:15, pN = 0.25 , pK = top_pK, nExp = calc_nExp, reuse.pANN = FALSE, sct = TRUE)


md <- sobj@meta.data
df_col_nm <-  grep("DF.classifications",colnames(md), value = T)
cell_num <- length(colnames(sobj)) #number of barcodes
dbs_plt <- DimPlot(sobj,pt.size = 1,label=TRUE, label.size = 5,reduction = "umap",group.by = df_col_nm )+theme(aspect.ratio = 1) + ggtitle(paste0(df_col_nm , "\n", "n=", cell_num))


print("Creating doublet plots...")
dir.create(path= paste0(args$out,"doublet_plots"), showWarnings= FALSE)
pdf(file=paste0(args$out,"doublet_plots/",args$lib_rna_atac_id, "_dbplots.pdf"), title= paste0(args$lib_rna_atac_id,"_doublet_plots"), width= 12, height= 10, compress=FALSE)
dbs_plt
dev.off()
#saveRDS(sobj_seurat_DF, "/projects/ps-epigen/users/kdang/eye_final/rds_objs/jb185_seurat_DF_0.05.rds")


}

## Add metadata columns

# library_id

# 2
#sobj[["library_rna_id"]] <- lib_rna_id
# adams #
#sobj[["sample_id"]] <- samp_id

#sobj[["tissue"]] <- tissue

#sobj[["cage"]] <- cage

#sobj[["sex"]] <- sex

#sobj[["age"]] <- age

#sobj[["treatment"]] <- treatment

#sobj[["treatment_a_b"]] <- treatment_a_b

#sobj[["age_days"]] <- age_days

# 11
#sobj[["birthdate"]] <- birthdate

#sobj[["library_atac_id"]] <- lib_atac_id

#sobj[["library_rna_atac_id"]] <- lib_rna_atac_id
# sobj[[deparse(substitute(doo))]] <- doo # make the name of the column the same as the variable name (doo)

# Saving RDS
dir.create(path= paste0(args$out,"df_objs/"), showWarnings= FALSE)
saveRDS(sobj, file = paste0(args$out, "df_objs/",args$lib_rna_atac_id,"_DF.0.05.RDS"))

# Saving RDS
#dir.create(path= paste0(args$out,"filtered_objs/"), showWarnings= FALSE)
#saveRDS(sobj, file = paste0(args$out, "filtered_objs/",lib_rna_atac_id,"_filt.RDS"))
