library("devtools")
load_all("~/artdeco/")

#library("artdeco")

inputfile = "~/MAPTor_NET/BAMs/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv"
visualization_data_path = "~/Deko/Data/Groetzinger_Scarpa.TPM.filtered.HGNC.Voom.TMM.normalized.tsv"
transcriptome_file_path = visualization_data_path
#visualization_data_path = system.file(
#    "/Data/Expression_data/Visualization_PANnen.tsv",
#    package = "artdeco"
#)

deconvolution_models = c(
    #"Four_differentiation_stages_Lawlor",
    #"Four_differentiation_stages_Lawlor_Progenitor",
    "Four_differentiation_stages_Segerstolpe_Prog_Hisc"
)

deconvolution_results_relative = artdeco::Determine_differentiation_stage(
    transcriptome_file_path = inputfile,
    models = "Alpha_Beta_Gamma_Delta_Baron_Progenitor_Stanescu_HESC_Yan"    
)

create_heatmap_differentiation_stages(
    transcriptome_file_path = visualization_data_path,
    deconvolution_results = deconvolution_results_relative,
    annotation_columns = c(
        "Differentiation_Stages_Subtypes",
        "Differentiation_Stages_Aggregated",
        "Differentiatedness"
    ),
    Graphics_parameters = "",
    baseline = "relative"
)

transcriptome_file_path = visualization_data_path
deconvolution_results = deconvolution_results_relative
annotation_columns = c(
    "Differentiation_Stages_Subtypes",
    "Differentiation_Stages_Aggregated",
    "Differentiatedness"
)
Graphics_parameters = ""
baseline = "relative"
