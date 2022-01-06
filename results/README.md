# Analysis

- Rmd: Interactive analysis  
- README: brief description of each file

ambient_rna_analysis.Rmd - Ambient rna estimations to use to filter contamination in DE analysis http://bioconductor.org/books/3.13/OSCA.multisample/ambient-problems.html.
cell_cycle_analysis.Rmd - Cell cycle analysis using cyclone http://bioconductor.org/books/3.14/OSCA.advanced/cell-cycle-assignment.html
check_litter1_no_integration.Rmd - UMAP, clustering and markers for only litter1 without integration between WT and HARE5  
check_new_runs_11nov21.Rmd - Results of downsampling seurat objects to see if the average RNA levels per cell drives clustering 
check_new_runs_8nov21.Rmd - Check runs obtained with automated seurat pipeline
cluster_level_comparison_11nov21.Rmd - Results of downsampling seurat objects to see if the average RNA levels per cell drives clustering, at the clustering level
count_matrix_qc.Rmd - Notes on interesting papers that do a similar experiment
de_pseudobulk_analysis.Rmd - Pseudobulk analysis using scran and edgeR http://bioconductor.org/books/3.14/OSCA.multisample/multi-sample-comparisons.html
diff_abundance_analysis.Rmd - Detect changes in cell type composition http://bioconductor.org/books/3.14/OSCA.multisample/differential-abundance.html 
differential_expression_analysis.Rmd - Exploration of differentially expressed genes following Seurat
effect_of_malat_on_clusters.Rmd - Malat1 correlates with low cell quality, litter1 specific clusters are have Malat1
extreme_pseudobulk_analysis.Rmd - Pseudobulk analysis collapsing at the level of samples instead of clusters 
extreme_pseudobulk_batch.Rmd - Pseudobulk analysis collapsing at the level of samples instead of clusters but comparing between litters instead of genotype. I tried this to see if the DEGs were useful for qc
feature_plots_of_markers.Rmd - Feature plots of conserved markers
find_conserved_markers.Rmd - Find markers using integrated dataset
intermediate_progenitors_analysis.Rmd - Filter cells by tbr2 presence and re-do clustering , using only litter2
intermediate_progenitors_both_litters_analysis.Rmd - Filter cells by tbr2 presence and re-do clustering but with both litters to be able to do Differential Expression
litter2_analyses.Rmd - Re-do all analyses using only litter 2
litter2_tbr2_independent_clusters.Rmd - Cluster tbr2+ cells from WT and HARE5 independently
milo_analysis.Rmd - Differential Abundance using an alternative approach called MILO https://www.nature.com/articles/s41587-021-01033-z
milo_analysis_no_integration.Rmd - MILO without integration 
qc_of_batch1_clusters.Rmd - Check markers for cell death and stress to QC litter 1 specific clusters
re_do_clusters.Rmd - Re-do clusters from full dataset using the integrated slot instead of SCT
separate_clustering.Rmd - Cluster cells from WT and HARE5 independently
subcluster_tbr2.Rmd - Cluster all cells from WT and HARE5 independently
wnt_pathway_exploration.Rmd - Check if genes from the WNT pathway are differentially expressed between genotypes


# Tutorials
milo_tutorial.Rmd
ambient_rna_tutorial.Rmd
cell_cycle_tutorial.Rmd
de_pseudobulk.Rmd 
diff_abundance.Rmd

# Exploration of Alejo's and AJ's code
scMm_analysis_aj_with_mito.Rmd
umap_aj_pipeline.Rmd
umap_aj_pipeline_with_mito.Rmd
umap_alejo_pipeline.Rmd
umap_alejo_pipeline_with_filters.Rmd
