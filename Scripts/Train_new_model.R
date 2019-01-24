library(bseqsc)

expression_training_mat = read.table(
    "~/Deko/Data/Alpha_Beta_Gamma_Delta_Baron_Progenitor.tsv",
    sep ="\t",
    header = T,
    stringsAsFactors = F
)
colnames(expression_training_mat) = 
    str_replace(colnames(expression_training_mat), pattern = "\\.", "_")
colnames(expression_training_mat) =
    str_replace(colnames(expression_training_mat), pattern = "^X", "")

meta_data = meta_info[colnames(expression_training_mat),]
subtype_vector = str_to_lower(meta_data$Subtype)
Marker_Gene_List = list()

### Data cleansing

row_var = apply(expression_training_mat, FUN = var, MARGIN = 1)
col_var = apply(expression_training_mat, FUN = var, MARGIN = 2)
expression_training_mat = expression_training_mat[row_var != 0,col_var != 0]
expression_training_mat = expression_training_mat[rowSums(expression_training_mat) >= 1,]

training_nr_marker_genes = 100

source("~/Deko/Scripts/Marker_gene_identification.R")
for( subtype in unique(subtype_vector) ){
    Marker_Gene_List[[subtype]] = identify_marker_genes(
        expression_training_mat,
        subtype_vector,
        subtype,
        training_nr_marker_genes
    )
}
print("Finished extracting marker genes for subtypes")

# Prepare bseq training

