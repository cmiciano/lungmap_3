library(Seurat)
library(Signac)
library(parallel)
library(GenomicRanges)

read_objects <- function(ipath1, ipath2=NULL){
    objs <- list.files(ipath1)
    obj_list <- c()

    for(o in objs){
        temp_obj <- readRDS(paste0(ipath1,o))
        obj_list <- append(obj_list, temp_obj)
        # since adding metadata for lung removed barcodes from rownames, readd so createchromassay works for cells param
        rownames(temp_obj@meta.data) <- temp_obj$BARCODE
    }

    if (is.null(ipath2) == FALSE){
        objs2 <- list.files(ipath2)
        for(o in objs2){
        temp_obj2 <- readRDS(paste0(ipath2,o))

        obj_list <- append(obj_list, temp_obj2)
        }
    }

    return(obj_list)
}

# Prepare objects for merging
merge_prep <- function(sobj){
    #### Extract fragments as frags objects from each dataset
    frags <- Fragments(sobj)
    #### Built feature matrices
    counts <- FeatureMatrix(
      fragments = frags,
      features = combined.peaks,
      cells = rownames(sobj@meta.data)
    )
    #### Create Chromatin assays.. obviously
    assay <- CreateChromatinAssay(counts, fragments = frags, min.features = -1 )
    #### Create new object from each assay
    # ******* change project to orig.ident if appropriate *******
    sobj[["ATAC_comb"]] <- assay
    return(sobj)
}

# NOTE may need to change
extract_samp_nm <- function(obs_list){
    samp_nms <- c()
    for (o in obs_list) {
        samp_nm <- as.character(o@meta.data$orig.ident[1])
        samp_nms <- append(samp_nms, samp_nm)
    }
    return(samp_nms)
}


args = commandArgs(trailingOnly=TRUE)
ipath1 <- args[1]
outp <- args[2]
proj <- args[3]
# min.features = -1 important 
# for cells with 0 features, end up in infinite loop without this if the cell numbers are different between rna and atac
#ipath1 <- "/projects/ps-epigen/users/cmiciano/SenNet_Multiome/07_DoubletFinder_Rep2/df_objs_md_liver/"


li <- read_objects(ipath1=ipath1)

# change default assay to atac
# assign to default cellranger ATAC peaks for each object, so the peaks are actually used when calling featurematrix
for (obj in 1:length(li)) {
   print(li[[obj]])
   DefaultAssay(li[[obj]]) <- 'ATAC'

}

for (obj in 1:length(li)) {
   print(li[[obj]])
}



### Collect metadata ******* cannot do when meta colnames or # of meta columns are different *********
#for (o in objects){
#  tmeta <- as.data.frame(o@meta.data)
#  metadata <- rbind(metadata, tmeta)
#}


### create list of genomic ranges. Change to "objects"
granges <- mclapply(mc.cores=detectCores()*.75, X= li, FUN=Signac::granges)


### concatonate them
grang.c <- c()


for (g in granges){grang.c <- append(grang.c, g)}
### Create a unified set of peaks to quantify in each dataset


combined.peaks <- reduce(grang.c)
### create proprocessed merged object list


Sys.time()
pmo_list <- lapply( X= li, FUN=merge_prep)
Sys.time()

print(pmo_list) # now they all have the same number of features, using the same combined.peaks


experiments <- extract_samp_nm(pmo_list)

#proj <- "230821_multiome_merge_liver"

print("Removing ATAC assay")
# remove atac assay
for (i in 1:length(pmo_list)) {
    print(i)
    DefaultAssay(pmo_list[[i]]) <- 'ATAC_comb'
    print(pmo_list[[i]])
    pmo_list[[i]][['ATAC']]<- NULL
}

print("Merging experiments")
# adding experiments (library_id) to keep cell names unique 
Sys.time()
all.objs <- merge(
  x = pmo_list[[1]],
  y = pmo_list[2:length(pmo_list)],
       add.cell.ids = experiments,
                project = proj
)

Sys.time()
#outp <- "/projects/ps-epigen/users/cmiciano/SenNet_Multiome/08_multiome_merge/"
print("Saving objects")
dir.create(outp)
outp_obj <- paste0(outp, proj, ".RDS")
#saveRDS(all.objs, "/projects/ps-epigen/users/cmiciano/SenNet_Multiome/08_multiome_merge/230821_multiome_merge_liver.rds"
saveRDS(all.objs, outp_obj)

