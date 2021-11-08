library(Seurat)
script_args <- commandArgs(trailingOnly = TRUE)

obj_path <- script_args[1]
print(obj_path)

obj_name <- tools::file_path_sans_ext(obj_path)
obj_out <- paste0(obj_name, "_integrated.rds")
markers_out <- paste0(obj_name, "_integrated_markers.rds")

Mm.merged <- readRDS(obj_path)

Mm.merged <- SCTransform(Mm.merged, variable.features.n = 5000, vars.to.regress = "percent.mt", verbose = FALSE)

Mm.merged <- RunPCA(object = Mm.merged, npcs = 300, features = VariableFeatures(object = Mm.merged))

Mm.merged <- FindNeighbors(object = Mm.merged, dims = 1:100)
Mm.merged <- FindClusters(object = Mm.merged, resolution = 1.1)
Mm.merged <- RunUMAP(Mm.merged, dims = 1:100)

ifnb.list <- SplitObject(Mm.merged, split.by = "sample")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 5000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE, npcs = 300)

immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:20)

# Seurat object output
Mm.markers <- FindAllMarkers(immune.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Markers output
saveRDS(immune.combined.sct, file = obj_out)
saveRDS(Mm.markers, file = markers_out)
