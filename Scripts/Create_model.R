library("devtools")
library("NMF")
load_all("~/artdeco/")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("stringr")
library("bseqsc")

expr_raw = read.table("~/Deko_Projekt/Data/Mouse_progenitor_pancreas_scRNA/senescence.tsv",sep ="\t", header = T, row.names = 1)
#expr_raw = readRDS("~/Downloads/Tosti.Seurat.normalized.S78048.RDS")
#expr_raw = read.table("~/Downloads/Feature_Barcode_rawCountMatrix_Filtered-YFP+_all-samples_QCed.csv",sep =",", header = T, row.names = 1)
#expr_raw = read.table("~/Downloads/Feature_Barcode_rawCountMatrix_Tuft+EEC_all-samples_QCed.csv",sep =",", header = T, row.names = 1)
colnames(expr_raw) = str_replace(colnames(expr_raw),"^X","")
expr_raw[1:5,1:5]

#meta_info = read.table("~/Deko_Projekt/Misc/Meta_information_scRNA.tsv", sep ="\t", header = T)
rownames(meta_info) = meta_info$Sample

#meta_info = read.table("~/Downloads/GSE172380_Cluster+CelltypeLabel_YFP+_all-samples_QCed.csv", sep =",", header = T)
#meta_info = read.table("~/Downloads/GSE172380_Cluster+CelltypeLabel_Tuft+EEC_all-samples_QCed.csv", sep =",", header = T)
meta_info = read.table("~/Deko_Projekt/Misc/Meta_information_scRNA.tsv", sep ="\t", header = T)
rownames(meta_info) = meta_info$Sample

table(meta_info$Cluster)
#meta_info$Sample = meta_info$X
#rownames(meta_info) = make.names(meta_info$Sample)
#meta_info$Cluster = meta_info$celltypeLabel

meta_data = meta_info[colnames(expr_raw),]
#meta_data$Cluster = make.names(meta_data$Cluster)
subtype_vector = meta_data$Cluster
table(subtype_vector)

#candidates = which(subtype_vector %in% c("alpha","beta","gamma","delta","acinar-s","acinar-reg+","acinar-i","ductal","muc5b+ ductal"))
#candidates = which(subtype_vector %in% c("Alpha","Beta","Gamma","Delta","Acinar-i","MUC5B+ Ductal","Acinar-REG+"))
#candidates = which(subtype_vector %in% c("Alpha","Beta","Gamma","Delta","Acinar","Ductal"))
#expr_raw = expr_raw[,candidates]
meta_data = meta_info[colnames(expr_raw),]
subtype_vector_reduced = meta_data$Cluster
table(subtype_vector_reduced)

subtype_vector = as.character(str_remove_all(colnames(expr_raw), pattern = "_[0-9]*"))
table(subtype_vector_reduced)

amount_genes = 400
amount_samples = 300
model_name = "Senesys_ADR_ADROHT_PS"

selected_samples = c()

for ( cell_type in unique(subtype_vector)){
    coords = which(subtype_vector == cell_type )
    
    if (length(coords) >= amount_samples)
        coords = sample(coords, size = amount_samples)
    
    selected_samples = c(selected_samples, coords)
}
length(selected_samples)

expr = expr_raw[,selected_samples]
dim(expr)
#meta_data_reduced = meta_info[colnames(expr),]
subtype_vector_reduced = subtype_vector[selected_samples]
table(subtype_vector_reduced)

rownames(expr) = str_to_upper(rownames(expr))

add_deconvolution_training_model_bseqsc(
    transcriptome_data = expr, 
    model_name = model_name,
    subtype_vector =  subtype_vector_reduced,
    training_p_value_threshold = 0.05,
    training_nr_permutations = 1000,
    training_nr_marker_genes = amount_genes
)
