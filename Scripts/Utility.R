
prepare_result_matrix_freestyle = function(
    meta_data,
    scale_values = TRUE,
    p_value_threshold = 0.05
){
    cand_labels = c("Alpha","Beta","Gamma","Delta","Acinar","Ductal","HISC")
    for (cand in cand_labels){
        
        proportion = meta_data[cand]
    }

    return(meta_data)
}

#' remove_model
#'
#' \code{remove_model} Removes a trained model
#'
#' @param model_name Name of the model
#' @usage
#' remove_model(
#'     model_name
#' )
#' @examples
#' remove_model(
#'     model_name = "my_model",
#' )
#' @export
remove_model = function(
    model_name
){
    
    model_path = paste(model_name, ".RDS", sep = "")
    model_path = paste(
        system.file("Models/NMF", package = "artdeco"),
        model_path,
        sep = "/"
    )
    
    if (! file.exists(model_path) )
        stop(paste0("Did not find model with name ",model_name))
    
    if ( model_name != "my_model"){
        
        file.remove(model_path)
        print(paste0("Successfully removed model with name ",model_name))
    } else {
        print("Removed my_model successfully!")
    }
}
