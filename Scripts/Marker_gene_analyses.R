library("stringr")

schlesinger_t = read.table("~/Downloads/Schlesinger.tsv", header = T, row.names = 1, sep ="\t")
schlesinger_genes = rownames(schlesinger_t)

baron = readRDS("~/Downloads/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron.RDS")
marker_genes_baron = baron[[3]]
marker_genes_baron_ductal = marker_genes_baron$ductal

(sum(marker_genes_baron_ductal %in% schlesinger_genes == TRUE) / length(marker_genes_baron_ductal)) * 100

segerstolpe = readRDS("~/Downloads/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Segerstolpe.RDS")
marker_genes_segerstolpe = segerstolpe[[3]]
marker_genes_segerstolpe_ductal = marker_genes_segerstolpe$ductal

(sum(marker_genes_segerstolpe_ductal %in% schlesinger_genes == TRUE) / length(marker_genes_segerstolpe_ductal)) * 100

lawlor = readRDS("~/Downloads/Alpha_Beta_Gamma_Delta_Acinar_Ductal_Lawlor.RDS")
marker_genes_lawlor = lawlor[[3]]
marker_genes_lawlor_ductal = marker_genes_lawlor$ductal

(sum(marker_genes_lawlor_ductal %in% schlesinger_genes == TRUE) / length(marker_genes_lawlor_ductal)) * 100

tosti = readRDS("~/artdeco/inst/Models/bseqsc/Baron_Tosti_800_genes_100_samples_all_endocrine_all_exocrine_no_metaplastic.RDS")
marker_genes_tosti = tosti[[3]]
names(marker_genes_tosti)
marker_genes_tosti_ductal = marker_genes_tosti$ductal

(sum(marker_genes_tosti_ductal %in% schlesinger_genes == TRUE) / length(marker_genes_tosti_ductal)) * 100
length(marker_genes_tosti_ductal)

marker_genes_tosti_new = marker_genes_tosti$alpha
(sum(marker_genes_tosti_new %in% schlesinger_genes == TRUE) / length(marker_genes_tosti_new)) * 100

#### NEUROG 3

neurog = readRDS("~/Downloads/EEC_Neurog_3.RDS")
marker_genes_neurog = neurog[[3]]
names(marker_genes_neurog)
marker_genes_neurog = marker_genes_tosti$`EEC-Progenitor (Neurog3+)`

#write.table(marker_genes_tosti,"~/Downloads/Neurog3.tsv", sep ="\t", quote = F)
