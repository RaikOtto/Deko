
prepare_result_matrix = function(fits){
        for (nr_fit in 1:length(fits)){
            
            fit = fits[nr_fit]
            res_coeff = t(fit[[1]]$coefficients)
            res_cor   = fit[[1]]$stats
            
            res_coeff[ is.na(res_coeff) ] = 0.0
            res_cor[ is.na(res_cor) ] = 0.0
            
            p_value_threshold = 0.05
            not_sig_samples = rownames(res_cor)[res_cor[,"P-value"] > p_value_threshold]
            not_sig_samples
            
            if (nr_fit == 1){
            
                Alpha_sim_scalar            = log( res_coeff[,"alpha"] + 1 )
                meta_data$Alpha_sim_scalar  = rep("",length(Alpha_sim_scalar))
                meta_data$Alpha_sim         = rep("",length(Alpha_sim_scalar))
                meta_data$Alpha_sim_scalar  = Alpha_sim_scalar
                meta_data$Alpha_sim         = Alpha_sim_scalar
                meta_data$Alpha_sim[ Alpha_sim_scalar <= quantile(Alpha_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Alpha_sim[ Alpha_sim_scalar > quantile(Alpha_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Alpha_sim[ Alpha_sim_scalar > quantile(Alpha_sim_scalar, seq(0,1,.01)[67])] = "high"
                meta_data[not_sig_samples,"Alpha_sim"] = "not_sig"
                
                Beta_sim_scalar            = log( res_coeff[,"beta"] + 1 )
                meta_data$Beta_sim_scalar  = rep("",length(Beta_sim_scalar))
                meta_data$Beta_sim         = rep("",length(Beta_sim_scalar))
                meta_data$Beta_sim_scalar  = Beta_sim_scalar
                meta_data$Beta_sim         = Beta_sim_scalar
                meta_data$Beta_sim[ Beta_sim_scalar <= quantile(Beta_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Beta_sim[ Beta_sim_scalar > quantile(Beta_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Beta_sim[ Beta_sim_scalar > quantile(Beta_sim_scalar, seq(0,1,.01)[67])] = "high"
                meta_data[not_sig_samples,"Beta_sim"] = "not_sig"
                
                Gamma_sim_scalar            = log( res_coeff[,"gamma"] + 1 )
                meta_data$Gamma_sim_scalar  = rep("",length(Gamma_sim_scalar))
                meta_data$Gamma_sim         = rep("",length(Gamma_sim_scalar))
                meta_data$Gamma_sim_scalar  = Gamma_sim_scalar
                meta_data$Gamma_sim         = Gamma_sim_scalar
                meta_data$Gamma_sim[ Gamma_sim_scalar <= quantile(Gamma_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Gamma_sim[ Gamma_sim_scalar > quantile(Gamma_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Gamma_sim[ Gamma_sim_scalar > quantile(Gamma_sim_scalar, seq(0,1,.01)[67])] = "high"
                meta_data[not_sig_samples,"Gamma_sim"] = "not_sig"
                
                Delta_sim_scalar            = log( res_coeff[,"delta"] + 1 )
                meta_data$Delta_sim_scalar  = rep("",length(Delta_sim_scalar))
                meta_data$Delta_sim         = rep("",length(Delta_sim_scalar))
                meta_data$Delta_sim_scalar  = Delta_sim_scalar
                meta_data$Delta_sim         = Delta_sim_scalar
                meta_data$Delta_sim[ Delta_sim_scalar <= quantile(Delta_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Delta_sim[ Delta_sim_scalar > quantile(Delta_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Delta_sim[ Delta_sim_scalar > quantile(Delta_sim_scalar, seq(0,1,.01)[67])] = "high"
                meta_data[not_sig_samples,"Delta_sim"] = "not_sig"
                
                Differentiated_sim_scalar = Alpha_sim_scalar + Beta_sim_scalar + Gamma_sim_scalar + Delta_sim_scalar
                meta_data$Differentiated_sim_scalar     = rep("",length(Differentiated_sim_scalar))
                meta_data$Differentiated_sim            = rep("",length(Differentiated_sim_scalar))
                meta_data$Differentiated_sim_scalar     = Differentiated_sim_scalar
                meta_data$Differentiated_sim            = Differentiated_sim_scalar
                meta_data$Differentiated_sim[ Differentiated_sim_scalar <= quantile(Differentiated_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Differentiated_sim[ Differentiated_sim_scalar > quantile(Differentiated_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Differentiated_sim[ Differentiated_sim_scalar > quantile(Differentiated_sim_scalar, seq(0,1,.01)[67])] = "high"
                meta_data[not_sig_samples,"Differentiated_sim"] = "not_sig"
            
            }
            
            if ( nr_fit == 2 ){
            
                Prog_sim_scalar            = log( res_coeff[,"progenitor"] + 1 )
                meta_data$Prog_sim_scalar  = rep("",length(Prog_sim_scalar))
                meta_data$Prog_sim         = rep("",length(Prog_sim_scalar))
                meta_data$Prog_sim_scalar  = Prog_sim_scalar
                meta_data$Prog_sim         = Prog_sim_scalar
                meta_data$Prog_sim[ Prog_sim_scalar <= quantile(Prog_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Prog_sim[ Prog_sim_scalar > quantile(Prog_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Prog_sim[ Prog_sim_scalar > quantile(Prog_sim_scalar, seq(0,1,.01)[67])] = "high"
                meta_data[not_sig_samples,"Prog_sim"] = "not_sig"
            }
            
            if ( nr_fit == 3 ){
                
                Hisc_sim_scalar            = log( res_coeff[,"hisc"] + 1 )
                meta_data$Hisc_sim_scalar  = rep("",length(Hisc_sim_scalar))
                meta_data$Hisc_sim         = rep("",length(Hisc_sim_scalar))
                meta_data$Hisc_sim_scalar  = Hisc_sim_scalar
                meta_data$Hisc_sim         = Hisc_sim_scalar
                meta_data$Hisc_sim[ Hisc_sim_scalar <= quantile(Hisc_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Hisc_sim[ Hisc_sim_scalar > quantile(Hisc_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Hisc_sim[ Hisc_sim_scalar > quantile(Hisc_sim_scalar, seq(0,1,.01)[67])] = "high"
                meta_data[not_sig_samples,"Hisc_sim"] = "not_sig"
                
                maxi <<- apply(  res_coeff, FUN = which.max, MARGIN = 1 )
                meta_data$Differentiation_type = colnames(res_coeff)[maxi]
                meta_data$Differentiation_type[ meta_data$Differentiation_type %in% c("alpha","Alpha","beta","Beta","gamma","Gamma","delta","Gamma")] = "Differentiated"
                meta_data$Differentiation_type[ meta_data$Differentiation_type %in% c("progenitor")] = "Progenitor"
                meta_data$Differentiation_type[ meta_data$Differentiation_type %in% c("hisc")] = "HISC"
            }
            
            if ( nr_fit == 4 ){
                
                Acinar_sim_scalar = log( res_coeff[,"acinar"] + 1 )
                meta_data$Acinar_sim_scalar     = rep("",length(Acinar_sim_scalar))
                meta_data$Acinar_sim            = rep("",length(Acinar_sim_scalar))
                meta_data$Acinar_sim_scalar     = Acinar_sim_scalar
                meta_data$Acinar_sim            = Acinar_sim_scalar
                meta_data$Acinar_sim[ Acinar_sim_scalar <= quantile(Acinar_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Acinar_sim[ Acinar_sim_scalar > quantile(Acinar_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Acinar_sim[ Acinar_sim_scalar > quantile(Acinar_sim_scalar, seq(0,1,.01)[67])] = "high"
                meta_data[not_sig_samples,"Acinar_sim"] = "not_sig"
                
                Ductal_sim_scalar = log( res_coeff[,"ductal"] + 1 )
                meta_data$Ductal_sim_scalar     = rep("",length(Ductal_sim_scalar))
                meta_data$Ductal_sim            = rep("",length(Ductal_sim_scalar))
                meta_data$Ductal_sim_scalar     = Ductal_sim_scalar
                meta_data$Ductal_sim            = Ductal_sim_scalar
                meta_data$Ductal_sim[ Ductal_sim_scalar <= quantile(Ductal_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Ductal_sim[ Ductal_sim_scalar > quantile(Ductal_sim_scalar, seq(0,1,.01)[34])] = "low"
                meta_data$Ductal_sim[ Ductal_sim_scalar > quantile(Ductal_sim_scalar, seq(0,1,.01)[67])] = "high"
                meta_data[not_sig_samples,"Ductal_sim"] = "not_sig"
            }
        }

#max_mat = meta_data[c("Hisc_sim_scalar","Prog_sim_scalar","Differentiated_sim_scalar")]
#maxi = apply(  max_mat, FUN = which.max, MARGIN = 1 )
#meta_data$Differentiation_type = c("Hisc","Progenitor","Differentiated")[maxi]
#not_sig_samples = p_value_mat[meta_data$Name,]
#for(i in 1:length(meta_data$Differentiation_type)){
#    if( p_value_mat[i,maxi[i]] == "not_sig" )
#        meta_data[i,"Differentiation_type"] = "not_sig"
#}
}