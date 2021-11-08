library(ggplot2)
library(tiff)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
library(Signac)
options(Seurat.memsafe = TRUE)
# Read in data
data <- Read10X("data_raw/run_outs/run_count_G_A1_outs/filtered_feature_bc_matrix", gene.column =1) # After a lot of work, I realized the gene column is critical to load the file into R

wt1 <- CreateSeuratObject(
  counts=data,
  project="wt1",
  min.cells = 3,
  min.features = 500
)



data <- Read10X("data_raw/run_outs/run_count_G_B1_outs/filtered_feature_bc_matrix", gene.column =1)

ki1 <- CreateSeuratObject(
  counts=data,
  project="ki1",
  min.cells = 3,
  min.features= 500
)


data <- Read10X("data_raw/run_outs/run_count_G_C2_outs/filtered_feature_bc_matrix", gene.column =1)

wt2 <- CreateSeuratObject(
  counts=data,
  project="wt2",
  min.cells = 3,
  min.features=500
) 

data <- Read10X("data_raw/run_outs/run_count_G_D2_outs/filtered_feature_bc_matrix", gene.column =1)

ki2 <- CreateSeuratObject(
  counts=data,
  project="ki2",
  min.cells = 3,
  min.features=500
) 


# Add groups

wt1$group <- "WT"
wt2$group <- "WT"
ki1$group <- "KI"
ki2$group <- "KI"


# Merge
wt_vs_ki <- merge(wt1, y= c(  wt2,ki1,ki2),
                            add.cell.ids = c("WT", "WT",
                                             "KI", "KI"),
                            project = "wt_vs_ki")
wt1 <- NULL
wt2 <- NULL
ki1 <- NULL
ki2 <- NULL
data <- NULL
gc()

wt_vs_ki$group <- factor(wt_vs_ki$group,  levels = c('WT', 'KI'))

#Normalizing the data
wt_vs_ki <- NormalizeData(wt_vs_ki)

#Identification of highly variable features (feature selection)
wt_vs_ki <- FindVariableFeatures(wt_vs_ki, selection.method = "vst", nfeatures = 2000)

#Scaling the data

wt_vs_ki <- ScaleData(wt_vs_ki)


#Perform linear dimensional reduction
wt_vs_ki <- RunPCA(wt_vs_ki, features = VariableFeatures(object = wt_vs_ki))

saveRDS(wt_vs_ki,"data_output/surat_objects/wt_vs_ki_pca.rds")
wt_vs_ki <- readRDS("data_output/surat_objects/wt_vs_ki_pca.rds")

# Examine and visualize PCA results a few different ways
ElbowPlot(wt_vs_ki,ndims = 50) +
  geom_vline(xintercept = 40, color = "red")
  

#Cluster the cells


wt_vs_ki <- FindNeighbors(wt_vs_ki, dims = 1:40)
wt_vs_ki <- FindClusters(wt_vs_ki, resolution = 1.1)
#Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
wt_vs_ki <- RunUMAP(wt_vs_ki, dims = 1:30)

DimPlot(wt_vs_ki, reduction = "umap")
DimPlot(wt_vs_ki, reduction = "umap", group.by = "group")


#Finding differentially expressed features (cluster biomarkers)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones

wt_vs_ki_markers <- FindAllMarkers(wt_vs_ki, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(wt_vs_ki_markers, file = "data_output/surat_objects/wt_vs_ki_markers.rds")


VlnPlot(wt_vs_ki, features = c("ENSMUSG00000036904"))

FeaturePlot(wt_vs_ki, features = c("ENSMUSG00000036904"),
            label = T,
            label.size = 2)




