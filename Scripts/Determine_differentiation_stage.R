library("devtools")
load_all("~/artdeco/")

#library("artdeco")

transcriptome_file_path = inputfile = "~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv"
visualization_data_path = "~/Deko/Data/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv"


#visualization_data_path = "~/Deko/Data/Visualization_PANnen.tsv"
#transcriptome_file_path = "~/Deko/Data/PANnen_Test_Data.tsv"

transcriptome_mat_vis = read.table(
    visualization_data_path,
    sep="\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    row.names = 1
)
colnames(transcriptome_mat_vis) = str_replace_all(
    colnames(transcriptome_mat_vis),
    pattern = "^X",
    ""
)

not_sig_samples = rownames(deconvolution_results)[
    deconvolution_results$P_value >= p_value]
for( sample_name in not_sig_samples)
    deconvolution_results[
        sample_name,
        grep(
            colnames(deconvolution_results),
            pattern = 'percent|Sig_score|log_odds|P_value|model',
            invert = TRUE
        )
        ] = "Not_significant"

correlation_matrix = cor(transcriptome_mat_vis)
pcr = prcomp(t(correlation_matrix))


deconvolution_results = Determine_differentiation_stage(
    transcriptome_file_path = transcriptome_file_path,
    models = c("Alpha_Beta_Gamma_Delta_Lawlor","Progenitor_Stanescu_HISC_Haber")
)

create_heatmap_differentiation_stages(
    visualization_data_path = visualization_data_path,
    deconvolution_results = deconvolution_results,
    aggregate_differentiated_stages = F,
    relative_baseline = T,
    Graphics_parameters = ""
)

pheatmap::pheatmap(
    correlation_matrix,
    annotation_col = deconvolution_results[,c("alpha_similarity_relative","beta_similarity_relative","gamma_similarity_relative","delta_similarity_relative")],
    annotation_colors = Graphics_parameters,
    annotation_legend = TRUE,
    treeheight_col = 0,
    treeheight_row = 0,
    show_colnames = FALSE,
    show_rownames = FALSE
)

### add models

scRNA_file_path = "~/Deko/Data/Progenitor_Stanescu_HESC_Yan.tsv"
model_name = str_replace_all(scRNA_file_path,pattern = "\\.tsv","")
model_name = tail(str_split(model_name,pattern = "/")[[1]],1)

t = read.table(scRNA_file_path, sep ="\t", header = T, row.names = 1, nrows = 1)
colnames(t) = str_replace_all(colnames(t),pattern ="\\.","_")
colnames(t) = str_replace_all(colnames(t),pattern ="^X","")
meta_data = meta_info[colnames(t),]
subtype_vector = str_to_lower(meta_data$Subtype)
table(subtype_vector)

add_deconvolution_training_model(
    transcriptome_data_path = scRNA_file_path,
    model_name = model_name,
    str_to_lower(subtype_vector),
    training_p_value_threshold = 0.05,
    training_nr_permutations = 0,
    training_nr_marker_genes = 800
)
