
prepare_result_matrix_freestyle = function(
    meta_data,
    scale_values = TRUE,
    p_value_threshold = 0.05
){
        for (nr_fit in 1:2){
            
            if (nr_fit == 1){
                
                res_coeff = fits_coeff[[nr_fit]]
                res_cor   = fits_stats[[nr_fit]]
                colnames(res_coeff) = str_replace_all(colnames(res_coeff) ,"^X","")
                rownames(res_cor) = str_replace_all(rownames(res_cor) ,"^X","")
                
                res_coeff[ is.na(res_coeff) ] = 0.0
                res_cor[ is.na(res_cor) ] = 0.0

                res_coeff = t(res_coeff)
                colnames(res_coeff) = str_to_lower(colnames(res_coeff))
                
                meta_data$Differentation_score = rep(0,nrow(meta_data))

                if (scale_values) {
                    Alpha_sim_scalar= exp( res_coeff[,"alpha"] )
                    Beta_sim_scalar = exp( res_coeff[,"beta"] )
                    Gamma_sim_scalar= exp( res_coeff[,"gamma"] )
                    Delta_sim_scalar= exp( res_coeff[,"delta"] )
                    meta_data$Differentation_score = exp(res_cor[,4])
                } else {
                    Alpha_sim_scalar= ( res_coeff[,"alpha"]  )
                    Beta_sim_scalar = ( res_coeff[,"beta"]  )
                    Gamma_sim_scalar= ( res_coeff[,"gamma"]  )
                    Delta_sim_scalar= ( res_coeff[,"delta"]  )
                    meta_data$Differentation_score = res_cor[,4]
                }
                
                meta_data[,"alpha"] = Alpha_sim_scalar
                meta_data[,"beta"] = Beta_sim_scalar
                meta_data[,"gamma"] = Gamma_sim_scalar
                meta_data[,"delta"] = Delta_sim_scalar
                    
                meta_data$Alpha_similarity= rep("low",nrow(meta_data))
                meta_data$Alpha_similarity[
                    Alpha_sim_scalar > quantile(Alpha_sim_scalar,
                    seq(0,1,.01)[67])
                ] = "high"

                meta_data$Beta_similarity= rep("low",nrow(meta_data))
                meta_data$Beta_similarity[ 
                    Beta_sim_scalar > quantile(Beta_sim_scalar,
                    seq(0,1,.01)[67])
                ] = "high"
                
                
                meta_data$Gamma_similarity= rep("low",nrow(meta_data))
                meta_data$Gamma_similarity[ 
                    Gamma_sim_scalar > quantile(Gamma_sim_scalar,
                    seq(0,1,.01)[67])
                ] = "high"
                
                meta_data$Delta_similarity= rep("low",nrow(meta_data))
                meta_data$Delta_similarity[ 
                    Delta_sim_scalar > quantile(Delta_sim_scalar,
                                                seq(0,1,.01)[67])
                ] = "high"
                
                # p value setting
                
                if ("acinar" %in% colnames(res_coeff)){
                    
                    if (log_values){
                        acinar_sim_scalar = log( res_coeff[,"acinar"] + 1 )
                        Ductal_sim_scalar  = log( res_coeff[,"ductal"] + 1 )
                    } else {
                        acinar_sim_scalar = res_coeff[,"acinar"]
                        Ductal_sim_scalar  = res_coeff[,"ductal"]
                    }
                    
                    
                    meta_data$acinar_similarity= rep("low",nrow(meta_data))
                    meta_data$acinar_similarity[ 
                        acinar_sim_scalar > quantile(acinar_sim_scalar,
                        seq(0,1,.01)[67])
                    ] = "high"
                    
                    meta_data$Ductal_similarity= rep("low",nrow(meta_data))
                    meta_data$Ductal_similarity[ 
                        Ductal_sim_scalar > quantile(Ductal_sim_scalar,
                        seq(0,1,.01)[67])
                    ] = "high"
                    
                    not_sig_index =
                        ( is.na(p_values) ) |
                        ( res_cor[sample,1] > p_value_threshold )
                        
                    meta_data[
                        not_sig_index,
                        grep(colnames(meta_data),
                        pattern = "acinar_sim|Ductal_sim")
                    ] = "not_sig"
                    }
                
                
                Differentiated_sim_scalar = 
                    Alpha_sim_scalar +
                    Beta_sim_scalar +
                    Gamma_sim_scalar +
                    Delta_sim_scalar
                meta_data$Differentiated_similarity = rep("low",nrow(meta_data))
                Differentiated_sim_scalar = exp(Differentiated_sim_scalar)
                meta_data$Differentiated_similarity[
                    Differentiated_sim_scalar > quantile(
                        Differentiated_sim_scalar,
                        seq(0,1,.01)[50]
                    )
                ] = "high"
                
                # set p-values
                
                p_values = res_cor[,1]
                p_values[is.na(p_values)]  = 1
                p_values[p_values == 9999] = 0
                no_sig_index = 
                    p_values >= p_value_threshold
                
                index = colnames(meta_data)[
                    str_to_lower(colnames(meta_data)) %in%
                        c(
                            "alpha_similarity",
                            "beta_similarity",
                            "gamma_similarity",
                            "delta_similarity",
                            "acinar_similarity",
                            "ductal_similarity"
                        )
                ]
                
                meta_data[no_sig_index,index] = "not_sig"
                }
            ## end fit 1
            
            if ( nr_fit == 2 ){
                
                res_coeff = t(fits_coeff[[nr_fit]])
                res_cor   = fits_stats[[nr_fit]]
                colnames(res_coeff) = str_replace_all(colnames(res_coeff) ,"^X","")
                rownames(res_cor) = str_replace_all(rownames(res_cor) ,"^X","")
                
                res_coeff[ is.na(res_coeff) ] = 0.0
                res_cor[ is.na(res_cor) ] = 0.0
                
                meta_data$De_differentation_score = rep(0, nrow(meta_data))
                meta_data$De_differentation_score = res_cor[,4]
                
                colnames(res_coeff) = str_to_lower(colnames(res_coeff))
                #if ( sum( colnames(res_coeff) %in% c("hisc","hesc")) > 0 ) {
                #    colnames(res_coeff)[which(colnames(res_coeff) %in% c("hisc","hesc"))] = "stem_cell"
                #}
                for ( label in colnames(res_coeff) ) {
                
                    similarity_label = paste(label,"similarity",sep="_")
                    if(scale_values){
                        sim_scalar = exp( res_coeff[,label] )
                        meta_data$De_differentation_score = exp(meta_data$De_differentation_score)
                    } else {
                        sim_scalar = ( res_coeff[,label] + 1 )
                    }
                    
                    meta_data[,label] = rep("low",nrow(meta_data))
                    meta_data[,label] = as.double(sim_scalar)
                    meta_data[,similarity_label]   = rep("low",nrow(meta_data))
                    meta_data[
                        sim_scalar > quantile(
                            sim_scalar,
                            seq(0,1,.01)[80]
                        ),similarity_label
                    ] = "high"
                }
                
                p_values = res_cor[,1]
                p_values[is.na(p_values)]  = 1
                p_values[p_values == 9999] = 0
                no_sig_index = 
                    p_values >= p_value_threshold
                
                index = colnames(meta_data)[
                    str_to_lower(colnames(meta_data)) %in%
                        c(
                            "hisc_similarity",
                            "hesc_similarity",
                            "progenitor_similarity"
                        )
                ]
                meta_data[no_sig_index,index] = "not_sig"
            }
            ## end fit 2
        }
    
    cand_labels_max = colnames(meta_data)[cand_label_indices]
    cand_labels_max = cand_labels_max[maxi]
    meta_data$Differentiation_stage = rep("",nrow(meta_data))
    meta_data$Differentiation_stage = cand_labels_max
    
    meta_data$Differentation_score[
        meta_data$Differentation_score >= quantile(meta_data$Differentation_score)[4]
    ] = quantile(meta_data$Differentation_score)[4]
    meta_data$Differentation_score[
        meta_data$Differentation_score <= quantile(meta_data$Differentation_score)[1]
        ] = quantile(meta_data$Differentation_score)[2]
    
    meta_data$De_differentation_score[
        meta_data$De_differentation_score >= quantile(meta_data$De_differentation_score)[4]
        ] = quantile(meta_data$De_differentation_score)[4]
    meta_data$De_differentation_score[
        meta_data$De_differentation_score <= quantile(meta_data$De_differentation_score)[1]
        ] = quantile(meta_data$De_differentation_score)[2]
    
    colnames(meta_data)[colnames(meta_data) == "progenitor_similarity"] = "Progenitor_similarity"
    colnames(meta_data)[colnames(meta_data) == "hisc_similarity"] = "HISC_similarity"
    colnames(meta_data)[colnames(meta_data) == "hesc_similarity"] = "HESC_similarity"
    
    colnames(meta_data)[colnames(meta_data) == "progenitor"] = "Progenitor"
    colnames(meta_data)[colnames(meta_data) == "hisc"] = "HISC"
    colnames(meta_data)[colnames(meta_data) == "hisc"] = "HESC"
    colnames(meta_data)[colnames(meta_data) == "alpha"] = "Alpha"
    colnames(meta_data)[colnames(meta_data) == "beta"] = "Beta"
    colnames(meta_data)[colnames(meta_data) == "gamma"] = "Gamma"
    colnames(meta_data)[colnames(meta_data) == "delta"] = "Delta"
    colnames(meta_data)[colnames(meta_data) == "acinar"] = "acinar"
    colnames(meta_data)[colnames(meta_data) == "ductal"] = "Ductal"
    
    cand_labels = c("alpha","beta","gamma","delta","hisc","hesc","progenitor")
    cand_label_indices = which( str_to_lower(colnames(meta_data)) %in% cand_labels)
    
    maxi = apply( 
        meta_data[,
                  cand_label_indices
                  ],
        FUN = which.max,
        MARGIN = 1
        
    )
    
    #p_value differentiation stage
    index = colnames(meta_data)[
        str_to_lower(colnames(meta_data)) %in%
            c(
                "Differentiation_stage",
                "Differentiated_similarity"
            )
        ]
    
    no_sig_index = apply(meta_data,MARGIN =1, FUN = function(vec){
        return("not_sig" %in% vec)
    })
    
    meta_data[no_sig_index,index] = "not_sig"

    return(meta_data)
}
