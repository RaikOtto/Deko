
prepare_result_matrix = function(
    meta_data,
    scale_values = TRUE
){
        for (nr_fit in 1:2){
            
            if (nr_fit == 1){
                
                res_coeff = fits_coeff[[nr_fit]]
                res_cor   = fits_stats[[nr_fit]]
                colnames(res_coeff) = str_replace_all(colnames(res_coeff) ,"^X","")
                rownames(res_cor) = str_replace_all(rownames(res_cor) ,"^X","")
                
                res_coeff[ is.na(res_coeff) ] = 0.0
                res_cor[ is.na(res_cor) ] = 0.0
                
                p_value_threshold = 0.5
                not_sig_samples = rownames(res_cor)[res_cor[,"P-value"] > p_value_threshold]
                not_sig_samples
                
                res_coeff = t(res_coeff)
                colnames(res_coeff) = str_to_lower(colnames(res_coeff))
            
                if (scale_values) {
                    Alpha_sim_scalar= exp( res_coeff[,"alpha"] )
                    Beta_sim_scalar = exp( res_coeff[,"beta"] )
                    Gamma_sim_scalar= exp( res_coeff[,"gamma"] )
                    Delta_sim_scalar= exp( res_coeff[,"delta"] )
                } else {
                    Alpha_sim_scalar= ( res_coeff[,"alpha"]  )
                    Beta_sim_scalar = ( res_coeff[,"beta"]  )
                    Gamma_sim_scalar= ( res_coeff[,"gamma"]  )
                    Delta_sim_scalar= ( res_coeff[,"delta"]  )
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
                
                if ("accinar" %in% colnames(res_coeff)){
                    
                    if (log_values){
                        Accinar_sim_scalar = log( res_coeff[,"accinar"] + 1 )
                        Ductal_sim_scalar  = log( res_coeff[,"ductal"] + 1 )
                    } else {
                        Accinar_sim_scalar = res_coeff[,"accinar"]
                        Ductal_sim_scalar  = res_coeff[,"ductal"]
                    }
                    
                    
                    meta_data$Accinar_similarity= rep("low",nrow(meta_data))
                    meta_data$Accinar_similarity[ 
                        Accinar_sim_scalar > quantile(Accinar_sim_scalar,
                        seq(0,1,.01)[67])
                    ] = "high"
                    
                    meta_data$Ductal_similarity= rep("low",nrow(meta_data))
                    meta_data$Ductal_similarity[ 
                        Ductal_sim_scalar > quantile(Ductal_sim_scalar,
                        seq(0,1,.01)[67])
                    ] = "high"
                }
                
                Differentiated_sim_scalar = Alpha_sim_scalar + Beta_sim_scalar + Gamma_sim_scalar + Delta_sim_scalar
                meta_data$Differentiated_similarity = rep("low",nrow(meta_data))
                Differentiated_sim_scalar = exp(Differentiated_sim_scalar)
                meta_data$Differentiated_similarity[
                    Differentiated_sim_scalar > quantile(
                        Differentiated_sim_scalar,
                        seq(0,1,.01)[50]
                    )
                ] = "high"
                
                meta_data$Differentation_score = rep(0,nrow(meta_data))
                meta_data$Differentation_score = res_cor[,4]
                meta_data$Differentation_score = log(meta_data$Differentation_score + 1)
                
                meta_data$P_value = rep(0,nrow(meta_data))
                meta_data$P_value = res_cor[,1]
            } ## end fit 1
            
            
            
            if ( nr_fit == 2 ){
                
                res_coeff = t(fits_coeff[[nr_fit]])
                res_cor   = fits_stats[[nr_fit]]
                colnames(res_coeff) = str_replace_all(colnames(res_coeff) ,"^X","")
                rownames(res_cor) = str_replace_all(rownames(res_cor) ,"^X","")
                
                res_coeff[ is.na(res_coeff) ] = 0.0
                res_cor[ is.na(res_cor) ] = 0.0
                
                p_value_threshold = 0.05
                not_sig_samples = rownames(res_cor)[res_cor[,"P-value"] > p_value_threshold]
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
                            seq(0,1,.01)[67]
                        ),similarity_label
                    ] = "high"
                }
            }
        }
    
    #cand_labels = c("alpha","beta","gamma","delta","hisc","hesc","progenitor")
    #cand_labels = cand_labels[cand_labels %in% colnames(meta_data)]
    #maxi = apply(  meta_data[,c("Differentation_score","De_differentation_score")], FUN = which.max, MARGIN = 1 )
    #meta_data$Differentiation_type = rep("",nrow(meta_data))
    #meta_data$Differentiation_type = colnames(meta_data[cand_labels])[maxi]
    #meta_data$Differentiation_type[ meta_data$Differentiation_type %in% c("alpha","Alpha","beta","Beta","gamma","Gamma","delta","Gamma")] = "Differentiated"
    #meta_data$Differentiation_type[ meta_data$Differentiation_type %in% c("progenitor")] = "Progenitor"
    #meta_data$Differentiation_type[ meta_data$Differentiation_type %in% c("hisc")] = "HISC"
    #meta_data$Differentiation_type[ meta_data$Differentiation_type %in% c("hesc")] = "HESC"
    
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
    colnames(meta_data)[colnames(meta_data) == "accinar"] = "Accinar"
    colnames(meta_data)[colnames(meta_data) == "ductal"] = "Ductal"
    
    meta_data$P_value = rep(0,nrow(meta_data))
    meta_data$P_value = res_cor[,1]

return(meta_data)
#max_mat = meta_data[c("Hisc_sim_scalar","Prog_sim_scalar","Differentiated_sim_scalar")]
#maxi = apply(  max_mat, FUN = which.max, MARGIN = 1 )
#meta_data$Differentiation_type = c("Hisc","Progenitor","Differentiated")[maxi]
#not_sig_samples = p_value_mat[meta_data$Name,]
#for(i in 1:length(meta_data$Differentiation_type)){
#    if( p_value_mat[i,maxi[i]] == "not_sig" )
#        meta_data[i,"Differentiation_type"] = "not_sig"
#}
}
