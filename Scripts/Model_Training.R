train_deconvolution_model = function(
    training_data = "",
    model_name = "",
    subtype_vector,
    training_p_value_threshold = 0.05,
    training_nr_permutations = 100,
    training_nr_marker_genes = 100
){
    
    library("stringr")
    library("bsseq")

    if( model_name == "")
        stop("Require model name, aborting")

    model_path = paste(c("~/ArtDeco/inst/Models/",model_name,".RDS"), collapse = "")
    
    if(file.exists(model_path))
        stop(paste0( collapse= "",
            c("Modelname ",model_name,
              " already exsits, please choose different name or delete existing model"))
        )

    if( ! file.exists(training_data)){
        stop(paste(
            c("Could not find file ",training_data,", aborting"),
            collapse = ""
        ))
    }

    expression_training_mat = read.table(
        training_data,
        sep ="\t",
        header = T,
        stringsAsFactors = F
    )
    colnames(expression_training_mat) = 
        str_replace(colnames(expression_training_mat), pattern = "\\.", "_")
    colnames(expression_training_mat) =
        str_replace(colnames(expression_training_mat), pattern = "^X", "")
    
    if (length(subtype_vector) == 0)
        stop(paste0("You have to provide the sample subtypes labels for model training"))
    subtype_vector = str_to_lower(subtype_vector)
    Marker_Gene_List = list()
    
    ### Data cleansing
    
    row_var = apply(expression_training_mat, FUN = var, MARGIN = 1)
    col_var = apply(expression_training_mat, FUN = var, MARGIN = 2)
    expression_training_mat = expression_training_mat[row_var != 0,col_var != 0]
    expression_training_mat = expression_training_mat[rowSums(expression_training_mat) >= 1,]
    
    source("~/Deko/Scripts/Dif_exp_script.R")
    for( subtype in unique(subtype_vector) ){
        Marker_Gene_List[[subtype]] = return_differentially_expressed_genes(
            expression_training_mat,
            subtype_vector,
            subtype,
            training_nr_marker_genes
        )
    }
    print("Finished extracting marker genes for subtypes")
    
    # Prepare bseq training
    
    training_mat_bseq = new(
        "ExpressionSet",
        exprs = as.matrix(expression_training_mat)
    )
    fData(training_mat_bseq) = data.frame( subtype_vector )
    pData(training_mat_bseq) = data.frame( subtype_vector )
    
    Basis = bseqsc_basis(
        training_mat_bseq,
        Marker_Gene_List,
        clusters = 'subtype_vector',
        samples = colnames(exprs(training_mat_bseq)),
        ct.scale = FALSE
    )
    
    print(
        "Basis trained, estimating deconvolution thresholds, this may take some time")
    
    test_mat = new(
        "ExpressionSet",
        exprs = as.matrix(expression_training_mat)
    );
    
    fit = bseqsc_proportions(
        test_mat,
        Basis,
        verbose = FALSE,
        absolute = T,
        log = F,
        perm = training_nr_permutations
    )
    
    print("Finished threshold determination")
    
    res_coeff = t(fit$coefficients)
    res_coeff_mat = as.double(unlist(res_coeff))
    res_coeff_mat = as.data.frame(
        matrix(
            res_coeff_mat,
            ncol = ncol(res_coeff),
            nrow = nrow(res_coeff)
        )
    )
    rownames(res_coeff_mat) = rownames(res_coeff)
    colnames(res_coeff_mat) = colnames(res_coeff)
    res_cor   = fit$stats
    
    res_coeff[ is.na(res_coeff) ] = 0.0
    res_cor[ is.na(res_cor) ] = 0.0
    
    self_scores = list()
    for (subtype in unique(subtype_vector)){
        self_scores[[subtype]] = as.double(
            res_coeff[
                which(subtype_vector == subtype),
                subtype
            ]
        )
    }
    
    model = list(Basis, self_scores, Marker_Gene_List)
    saveRDS(model,model_path)
    
    print(paste0("Storing model: ", model_path))
    
    print(paste0("Finished training model: ", model_name))
    
}
