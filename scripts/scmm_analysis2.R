library(Seurat)
script_args <- commandArgs(trailingOnly = TRUE)

obj_path <- script_args[1]
print(obj_path)

obj_name <- tools::file_path_sans_ext(obj_path)
obj_out <- paste0(obj_name, "_integrated.rds")
markers_out <- paste0(obj_name, "_integrated_markers.rds")

seurat_obj <- readRDS(obj_path)

seurat_obj_list <- SplitObject(seurat_obj, split.by = "sample")
seurat_obj_list <- lapply(X = seurat_obj_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = seurat_obj_list, nfeatures = 5000)
seurat_obj_list <- PrepSCTIntegration(object.list = seurat_obj_list, anchor.features = features)

seurat_obj <- FindIntegrationAnchors(object.list = seurat_obj_list, normalization.method = "SCT",
                                         anchor.features = features)
seurat_obj <- IntegrateData(anchorset = seurat_obj, normalization.method = "SCT")
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, npcs = 300)

seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)

seurat_obj <- FindNeighbors(object = seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(object = seurat_obj, resolution = 1.1)

# Seurat object output
seurat_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

# Markers output
saveRDS(seurat_obj, file = obj_out)
saveRDS(seurat_markers, file = markers_out)