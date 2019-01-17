source("~/Deko/Scripts/Utility.R")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

vis_data = read.table("~/Deko/Data/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",sep ="\t", header = T, stringsAsFactors = F,row.names = 1)
rownames(vis_data) = str_to_upper(rownames(vis_data))
colnames(vis_data) = str_replace_all(colnames(vis_data), "\\.", "_")
colnames(vis_data) = str_replace_all(colnames(vis_data), "^X", "")
cor_mat = cor(vis_data)

bam_data = read.table("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t", header = T, stringsAsFactors = F,row.names = 1)
rownames(bam_data) = str_to_upper(rownames(bam_data))
colnames(bam_data) = str_replace_all(colnames(bam_data), "\\.", "_")
colnames(bam_data) = str_replace_all(colnames(bam_data), "^X", "")

eset = new("ExpressionSet", exprs=as.matrix(bam_data[,]));

B_differentiated = readRDS("~/artdeco/inst/Models/Alpha_Beta_Gamma_Delta_Segerstolpe.RDS")[[1]]
B_undifferentiated = readRDS("~/artdeco/inst/Models/Progenitor_Stanescu_HESC_Yan_HISC_Haber.RDS")[[1]]


nr_permutations = 0
fit_differentiated = bseqsc_proportions(eset, B_differentiated, verbose = TRUE, absolute = T, perm = nr_permutations, log = F)
#fit_differentiated_six = bseqsc_proportions(eset, B_differentiated_six, verbose = FALSE, absolute = T, log = F, perm = nr_permutations)
fit_undifferentiated = bseqsc_proportions(eset, B_undifferentiated, verbose = TRUE, absolute = T, log = F, perm = nr_permutations)

meta_data = meta_info[colnames(exprs(eset)),c("NEC_NET","OS_Tissue")]
fits_coeff =list(fit_differentiated$coefficients,fit_undifferentiated$coefficients)
fits_stats =list(fit_differentiated$stats,fit_undifferentiated$stats)

meta_data = prepare_result_matrix(meta_data, scale_values = F)
#meta_data$Diff_Type =  c("Differentiated","Progenitor","HSC")[maxi]
vis_mat = meta_data[,c(
    "Differentiated_similarity",
    "Progenitor_similarity",
    "HISC_similarity",
    "HESC_similarity",
    "De_differentation_score",
    "Differentation_score"
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
meta_data[,c("De_differentation_score","Differentation_score")]=quantile_normalisation(meta_data[,c("De_differentation_score","Differentation_score")])
