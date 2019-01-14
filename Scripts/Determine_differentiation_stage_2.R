Determine_PANnen_differentiation_stage = function(
    inputfile = "~/MAPTor_NET/BAMs/TPMs.57_Samples.Groetzinger_Scarpa.Non_normalized.HGNC.tsv",
    deconvolution_models = c(
        "Four_differentiation_stages_Lawlor",
        "Four_differentiation_stages_Lawlor_Progenitor",
        "Four_differentiation_stages_Segerstolpe_Prog_Hisc"
    ),
    nr_permutations = 100
){
    
    if (!file.exists(inputfile)){
        stop(paste0("Could not find file ",inputfile))
    }
    
    inputfile = read.table(
        inputfile,
        sep="\t",
        header = T,
        stringsAsFactors = F
    )
    deconvolution_data = new("ExpressionSet", exprs=as.matrix(inputfile));

    if (length(deconvolution_models) == 0)
        stop("Require at least one model")
    
    model_list = list()
        
    for (model in deconvolution_models){
        model_path = paste0(c("~/ArtDeco/inst/Models/",model,".RDS"),collapse = "")
        if( !file.exists(model_path))
            stop(paste0(c("Could not find model ",model_path,", aborting"), collapse = ""))
        model_list[model] = readRDS(model_path)
    }
    print("Models loaded")
    
    prediction_list = list()
    
    for (model in deconvolution_models){
        print(past0("Deconvolving with model: ",model))
        model_instance = model_list[model]
        basis_mat = model
        fit = bseqsc_proportions(
            deconvolution_data,
            model,
            verbose = FALSE,
            absolute = T,
            log = F,
            perm = nr_permutations
        )
        prediction_list[model] = fit
    }
    
    fits = list(fit_differentiated, fit_progenitor, fit_undifferentiated)
    meta_data <<- meta_info[colnames(eset),]
    
    source("~/Deko/Scripts/Utility_script.R")
}
