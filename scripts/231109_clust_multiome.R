### rlibs
library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(parallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))


library(harmony)

# NOTE
sobj <- readRDS("/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231108_01_merged_multiome_lung8.RDS")

print(sobj)

print(table(sobj$orig.ident))

# Skipping SCT, doing log normalize instead
DefaultAssay(sobj) <- 'RNA'


pdf("/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231109_qc_nct_libid.pdf")
VlnPlot(sobj, "nCount_RNA", split.by = "orig.ident")
dev.off()


pdf("/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231109_qc_nft_libid.pdf")
VlnPlot(sobj, "nFeature_RNA", split.by = "orig.ident")
dev.off()

pdf("/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231109_qc_mito_libid.pdf")
VlnPlot(sobj, "percent.mt", split.by = "orig.ident")
dev.off()


Sys.time()
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj), verbose = F)
sobj <- FindNeighbors(sobj, dims = 1:10)
sobj <- FindClusters(sobj, resolution = 0.5)
sobj <- RunUMAP(sobj, dims = 1:10, verbose = F, reduction.name = "umap.rna")
Sys.time()


# NOTE
# non batch corrected umap
pdf("/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231109_umap_nobatch_rna_libid.pdf")
DimPlot(sobj, reduction = "umap.rna", group.by = "orig.ident")
dev.off()



#NOTE
pdf("/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231109_umap_nobatch_rna_clust.pdf")
DimPlot(sobj, reduction = "umap.rna", group.by = "seurat_clusters")
dev.off()



## ATAC processing
DefaultAssay(sobj) <- 'ATAC_comb'

### Processing
Sys.time()
sobj <- RunTFIDF(sobj) %>%
FindTopFeatures(min.cutoff = 'q0') %>%
RunSVD() 
Sys.time()

Sys.time()
sobj <- RunUMAP(sobj, reduction = 'lsi', dims = 2:30, verbose = F,
        reduction.name='umap.atac', reduction.key='atac_UMAP_')
Sys.time()



sobj <- FindNeighbors(sobj, reduction = 'lsi', dims = 2:30) %>%
FindClusters(algorithm = 3, resolution = .3, verbose = FALSE)



# NOTE
pdf("/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231109_umap_nobatch_atac_clust.pdf")
DimPlot(sobj, reduction = "umap.atac", group.by = "seurat_clusters")
dev.off()


# wnn using no batch reduction
Sys.time()
sobj <- FindMultiModalNeighbors(sobj, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:30),
            weighted.nn.name = "weighted.nn",
            snn.graph.name = "wsnn")
Sys.time()



Sys.time()
sobj <- RunUMAP(sobj, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnn_UMAP_',
verbose = F)
Sys.time()


# try algorithm 3 instead of 4 and method 'igraph'
sobj <- FindClusters(sobj, graph.name='wsnn', algorithm=3,  resolution = 0.25, verbose=FALSE)

# NOTE
pdf("/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231109_umap_nobatch_wnn_clust.pdf")
DimPlot(sobj, reduction = "umap.wnn", group.by = "seurat_clusters")
dev.off()


# plot seurat_clusters based on WNN 
options(repr.plot.width=18, repr.plot.height=6)
p1 <- DimPlot(sobj , reduction='umap.rna', group.by='seurat_clusters', label=TRUE, label.size=3, repel=TRUE) + ggtitle('RNA')
p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('RNA only')
p2 <- DimPlot(sobj , reduction='umap.atac', group.by='seurat_clusters', label=TRUE, label.size=3, repel=TRUE) + ggtitle('ATAC')
p2 <- p2 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('ATAC only')
p3 <- DimPlot(sobj , reduction='umap.wnn', group.by='seurat_clusters', label=TRUE, label.size=3, repel=TRUE) + ggtitle('WNN')
p3 <- p3 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('WNN - Combined')
p4 <- p1 + p2 + p3 & NoLegend() 

# NOTE
pdf("/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231109_umap_nobatch_all_modal.pdf", width = 15, height = 5)
print(p4)
dev.off()

# NOTE
Sys.time()
saveRDS(sobj, "/projects/ps-epigen/users/cmiciano/Lung/lungmap_3/05_multiome_md_merge_8/231109_02_clust_nobatch_raw.RDS")
Sys.time()


