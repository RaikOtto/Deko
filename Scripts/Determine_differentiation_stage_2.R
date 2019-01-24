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
#expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
expr_raw = read.table("~/artdeco/inst/Data/Expression_data/PANnen_Test_Data.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
expr = matrix(as.double(as.character(unlist(expr_raw[ rownames(expr_raw) %in% sad_genes,]))), ncol = ncol(expr_raw))
colnames(expr) = colnames(expr_raw)
rownames(expr) = rownames(expr_raw)[rownames(expr_raw) %in% sad_genes]
cor_mat = cor(expr[,]);pcr = prcomp(t(cor_mat))

bam_data = read.table("~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",sep ="\t", header = T, stringsAsFactors = F,row.names = 1)
rownames(bam_data) = str_to_upper(rownames(bam_data))
colnames(bam_data) = str_replace_all(colnames(bam_data), "\\.", "_")
colnames(bam_data) = str_replace_all(colnames(bam_data), "^X", "")


path_transcriptome_file = system.file(
    "/Data/Expression_data/PANnen_Test_Data.tsv",
    package="artdeco"
)
visualization_data_path = system.file(
    "/Data/Expression_data/Visualization_PANnen.tsv",
    package="artdeco")

path_transcriptome_file = "~/Deko/Data/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv"
visualization_data_path = "~/Deko/Data/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.subset.tsv"

#deconvolution_results_segerstolpe = deconvolution_results

deconvolution_results_baron = Determine_differentiation_stage(
    transcriptome_file_path = path_transcriptome_file,
    deconvolve_exokrine_tissue = FALSE,
    HISC_stem_cell_only = FALSE,
    models = c(
        "Alpha_Beta_Gamma_Delta_Baron",
        "Progenitor_Stanescu_HISC_Haber"
    ),
    nr_permutations = 100,
    output_file = ""
)

deconvolution_results_segerstolpe = Determine_differentiation_stage(
    transcriptome_file_path = path_transcriptome_file,
    deconvolve_exokrine_tissue = FALSE,
    HISC_stem_cell_only = FALSE,
    models = c(
        "Alpha_Beta_Gamma_Delta_Segerstolpe",
        "Progenitor_Stanescu_HESC_Yan_HISC_Haber"
    ),
    nr_permutations = 100,
    output_file = ""
)

deconvolution_results_lawlor = Determine_differentiation_stage(
    transcriptome_file_path = path_transcriptome_file,
    deconvolve_exokrine_tissue = FALSE,
    HISC_stem_cell_only = FALSE,
    models = c(
        "Alpha_Beta_Gamma_Delta_Lawlor",
        "Progenitor_Stanescu_HESC_Yan_HISC_Haber"
    ),
    nr_permutations = 100,
    output_file = ""
)



deconvolution_results_meta = deconvolution_results_segerstolpe

for(subtype in c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","Progenitor","HISC","HESC","Differentiation_score","De_differentiation_score","P_value_de_differentiated","P_value_differentiated")){
    if (!(subtype %in% colnames(deconvolution_results_meta)))
        next()
    deconvolution_results_meta[,subtype] =
        (deconvolution_results_segerstolpe[subtype] +  
        deconvolution_results_baron[subtype] + 
        deconvolution_results_lawlor[subtype])/3
}

deconvolution_results[,"NEC_NET"] = meta_info[rownames(deconvolution_results),"NEC_NET"]
deconvolution_results[,"Study"] = meta_info[rownames(deconvolution_results),"Study"]
vis_mat = create_heatmap_differentiation_stages(
    visualization_data_path,
    deconvolution_results,
    #high_threshold_diff = 99,
    #high_threshold_de_diff = 25,
    p_value_threshold = 33,
    aggregate_differentiated_stages = F
)
