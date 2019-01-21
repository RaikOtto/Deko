library(devtools)
load_all("~/artdeco")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

source("~/Deko/Scripts/Visualization_colors.R")
genes_of_interest_hgnc_t = read.table("~/Deko/Misc//Stem_signatures.gmt",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1

sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[13,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
#expr_raw = read.table("~/Deko/Data/Visualization_PANnen.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw))
colnames(expr) = colnames(expr_raw)
rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
cor_mat = cor(expr[,]);pcr = prcomp(t(cor_mat))

bam_data = read.table("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t", header = T, stringsAsFactors = F,row.names = 1)
rownames(bam_data) = str_to_upper(rownames(bam_data))
colnames(bam_data) = str_replace_all(colnames(bam_data), "\\.", "_")
colnames(bam_data) = str_replace_all(colnames(bam_data), "^X", "")

eset = new("ExpressionSet", exprs=as.matrix(bam_data[,]));

transcriptome_file_path = "~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv"
#transcriptome_file_path = "~/Deko/Data/TPMs.81_Samples.Groetzinger_Scarpa.tsv"

nr_permutations = 200
Determine_differentiation_stage(
    transcriptome_file_path
)

vis_mat = deconvolution_result[,c(
    "Differentiatedness",
    "Progenitor_similarity",
    "HISC_similarity",
    #"HESC_similarity",
    "Differentiation_score",
    "De_differentiation_score",
    "Correlation_de_differentiated"
    #"HESC_similarity",
)]
pheatmap::pheatmap(
    #t(res_coeff),
    cor_mat,
    #annotation_col = meta_data[c("Hisc_sim","Prog_sim","Differentiated_sim","Grading","NEC_NET")],
    annotation_col = vis_mat,
    annotation_colors = aka3,
    annotation_legend = T,
    treeheight_col = 0,
    treeheight_row = 0,
    show_colnames = F,
    show_rownames = F#,
    #color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(length(breaksList)),
    #cluster_cols = F, cluster_rows = F
)

quantile_normalisation <- function(df){
    df_rank <- apply(df,2,rank,ties.method="min")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)
    
    index_to_mean <- function(my_index, my_mean){
        return(my_mean[my_index])
    }
    
    df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
    rownames(df_final) <- rownames(df)
    return(df_final)
}

meta_data = deconvolution_result_relative
meta_data = meta_info[colnames((bam_data)),]
vis_mat = deconvolution_result_absolute
vis_mat = deconvolution_result_absolute[,c(
    "Alpha_similarity",
    "Beta_similarity",
    "Gamma_similarity",
    "Delta_similarity",
    "Differentiation_score"
)]
vis_mat$NEC_NET = meta_info[colnames(bam_data),"NEC_NET"]
vis_mat$Location = meta_info[colnames(bam_data),"Location"]
