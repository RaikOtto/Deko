library("devtools")
library("NMF")
load_all("~/artdeco/")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")

expr_raw = readRDS("~/Downloads/Tosti.Seurat.normalized.S78048.RDS")

meta_info = read.table("~/Deko_Projekt/Misc/Tosti_Metadaten.tsv", sep ="\t", header = T)
rownames(meta_info) = meta_info$Cell

meta_data = meta_info[colnames(expr_raw),]
subtype_vector = str_to_lower(meta_data$Cluster)
table(subtype_vector)

candidates = which(subtype_vector %in% c("alpha","beta","gamma","delta","acinar-s","acinar-reg+","acinar-i","ductal","muc5b+ ductal"))
#candidates = which(subtype_vector %in% c("alpha","beta","gamma","delta","acinar-s","ductal"))
expr_raw = expr_raw[,candidates]
meta_data = meta_info[colnames(expr_raw),]
subtype_vector_reduced = meta_data$Cluster
table(subtype_vector_reduced)

amount_genes = 300
amount_samples = 300
model_name = "Tosti_300_300_endocrine_exocrine_metaplastic"

selected_samples = c()

for ( cell_type in unique(subtype_vector_reduced)){
    coords = which(meta_data$Cluster == cell_type )
    
    if (length(coords) >= amount_samples)
        coords = sample(coords, size = amount_samples)
    
    selected_samples = c(selected_samples, coords)
}
length(selected_samples)

expr = expr_raw[,selected_samples]
dim(expr)
meta_data_reduced = meta_info[colnames(expr),]
table(meta_data_reduced$Cluster)

add_deconvolution_training_model_bseqsc(
    transcriptome_data = expr,
    model_name = model_name,
    subtype_vector =  str_to_lower(meta_data_reduced$Cluster),
    training_p_value_threshold = 0.05,
    training_nr_permutations = 0,
    training_nr_marker_genes = amount_genes
)

