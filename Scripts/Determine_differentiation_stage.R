library("devtools")
library("NMF")
load_all("~/artdeco/")

meta_info = read.table("~/Deco/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

### add models

scRNA_file_path = "~/Deco/Data/Alpha_Beta_Gamma_Delta_Hisc_Baron.tsv"
model_name = str_replace_all(scRNA_file_path,pattern = "\\.tsv","")
model_name = tail(str_split(model_name,pattern = "/")[[1]],1)

t = read.table(scRNA_file_path, sep ="\t", header = T, row.names = 1, nrows = 1)
colnames(t) = str_replace_all(colnames(t),pattern ="\\.","_")
colnames(t) = str_replace_all(colnames(t),pattern ="^X","")
dim(t)
meta_data = meta_info[colnames(t),]
subtype_vector = str_to_lower(meta_data$Subtype)
table(subtype_vector)

transcriptome_data = read.table(scRNA_file_path,sep="\t",header  = T)

add_deconvolution_training_model_bseqsc(
    transcriptome_data = transcriptome_data,
    model_name = model_name,
    subtype_vector =  str_to_lower(subtype_vector),
    training_p_value_threshold = 0.05,
    training_nr_permutations = 0,
    training_nr_marker_genes = 800
)

add_deconvolution_training_model_music(
    transcriptome_data = transcriptome_data,
    model_name = model_name,
    subtype_vector
)

add_deconvolution_training_model_NMF(
    transcriptome_data = transcriptome_data,
    model_name = model_name,
    subtype_vector = subtype_vector,
    parallel_processes = 8
)
