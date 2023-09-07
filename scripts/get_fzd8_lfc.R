# Filtering DEG results to get only Fzd8 lfc
all_clusters_deg_all <- readr::read_csv("data_output/surat_objects/all_clusters_deg_all.csv")
fzd8_clusters_deg_all <- readr::read_csv("data_output/surat_objects/fzd8_clusters_deg_all.csv")
tbr2_clusters_deg_all <- readr::read_csv("data_output/surat_objects/tbr2_clusters_deg_all.csv")


all_clusters_deg_all_filt <- all_clusters_deg_all %>%
  dplyr::filter(SYMBOL == "Fzd8")

fzd8_clusters_deg_all_filt <- fzd8_clusters_deg_all %>%
  dplyr::filter(SYMBOL == "Fzd8")

tbr2_clusters_deg_all_filt <- tbr2_clusters_deg_all %>%
  dplyr::filter(SYMBOL == "Fzd8")


readr::write_csv(all_clusters_deg_all_filt, "data_output/surat_objects/all_clusters_deg_all_fzd8.csv")
readr::write_csv(fzd8_clusters_deg_all_filt, "data_output/surat_objects/fzd8_clusters_deg_all_fzd8.csv")
readr::write_csv(tbr2_clusters_deg_all_filt, "data_output/surat_objects/tbr2_clusters_deg_all_fzd8.csv")