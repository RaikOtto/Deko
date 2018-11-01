
res_coeff = t(fit$coefficients)
res_cor   = fit$stats

res_coeff[ is.na(res_coeff) ] = 0.0
res_cor[ is.na(res_cor) ] = 0.0

meta_data = meta_info[ rownames(res_coeff),]

not_sig_samples = rownames(res_cor)[res_cor[,"P-value"] > .1]
not_sig_samples

alpha_sim = log( res_coeff[,"alpha"]+1)
meta_data$Alpha_sim = rep("",length(alpha_sim))
meta_data$Alpha_sim = alpha_sim
#meta_data$Alpha_sim[ alpha_sim <= quantile(alpha_sim, seq(0,1,.01)[34])] = "low"
#meta_data$Alpha_sim[ alpha_sim > quantile(alpha_sim, seq(0,1,.01)[34])] = "medium"
#meta_data$Alpha_sim[ alpha_sim > quantile(alpha_sim, seq(0,1,.01)[67])] = "high"
#meta_data[not_sig_samples,"Alpha_sim"] = "not_sig"

beta_sim = log( res_coeff[,"beta"]+1)
meta_data$Beta_sim = rep("",length(beta_sim))
meta_data$Beta_sim = beta_sim
#meta_data$Beta_sim[ beta_sim <= quantile(beta_sim, seq(0,1,.01)[34])] = "low"
#meta_data$Beta_sim[ beta_sim > quantile(beta_sim, seq(0,1,.01)[34])] = "medium"
#meta_data$Beta_sim[ beta_sim > quantile(beta_sim, seq(0,1,.01)[67])] = "high"
#meta_data[not_sig_samples,"Beta_sim"] = "not_sig"

gamma_sim = log( res_coeff[,"gamma"]+1)
meta_data$Gamma_sim = rep("",length(gamma_sim))
meta_data$Gamma_sim = gamma_sim
#meta_data$Gamma_sim[ gamma_sim <= quantile(gamma_sim, seq(0,1,.01)[34])] = "low"
#meta_data$Gamma_sim[ gamma_sim > quantile(gamma_sim, seq(0,1,.01)[34])] = "medium"
#meta_data$Gamma_sim[ gamma_sim > quantile(gamma_sim, seq(0,1,.01)[67])] = "high"
#meta_data[not_sig_samples,"Gamma_sim"] = "not_sig"

delta_sim = log( res_coeff[,"delta"]+1)
meta_data$Delta_sim = rep("",length(delta_sim))
meta_data$Delta_sim = delta_sim
#meta_data$Delta_sim[ delta_sim <= quantile(delta_sim, seq(0,1,.01)[34])] = "low"
#meta_data$Delta_sim[ delta_sim > quantile(delta_sim, seq(0,1,.01)[34])] = "medium"
#meta_data$Delta_sim[ delta_sim > quantile(delta_sim, seq(0,1,.01)[67])] = "high"
#meta_data[not_sig_samples,"Delta_sim"] = "not_sig"

acinar_sim = log( res_coeff[,"acinar"]+1)
meta_data$Acinar_sim = rep("",length(acinar_sim))
meta_data$Acinar_sim = acinar_sim
#meta_data$Acinar_sim[ acinar_sim <= quantile(acinar_sim, seq(0,1,.01)[34])] = "low"
#meta_data$Acinar_sim[ acinar_sim > quantile(acinar_sim, seq(0,1,.01)[34])] = "medium"
#meta_data$Acinar_sim[ acinar_sim > quantile(acinar_sim, seq(0,1,.01)[67])] = "high"
#meta_data[not_sig_samples,"Acinar_sim"] = "not_sig"

ductal_sim = log( res_coeff[,"ductal"]+1)
meta_data$Ductal_sim = rep("",length(ductal_sim))
meta_data$Ductal_sim = ductal_sim
#meta_data$Ductal_sim[ ductal_sim <= quantile(ductal_sim, seq(0,1,.01)[34])] = "low"
#meta_data$Ductal_sim[ ductal_sim > quantile(ductal_sim, seq(0,1,.01)[34])] = "medium"
#meta_data$Ductal_sim[ ductal_sim > quantile(ductal_sim, seq(0,1,.01)[67])] = "high"
#meta_data[not_sig_samples,"Ductal_sim"] = "not_sig"


diff_sim = log(rowSums(res_coeff[,c("alpha","beta","gamma","delta")])+1)
meta_data$Differentiated_sim = rep("low", length(diff_sim))
meta_data$Differentiated_sim[ diff_sim > quantile(diff_sim, seq(0,1,.01)[34])] = "medium"
meta_data$Differentiated_sim[ diff_sim > quantile(diff_sim, seq(0,1,.01)[67])] = "high"
meta_data[not_sig_samples,"Differentiated_sim"] = "not_sig"

prog_sim = log(res_coeff[,"e13.5"]+1)
meta_data$Progenitor_sim = rep("low", length(prog_sim))
meta_data$Progenitor_sim[ prog_sim > quantile(prog_sim, seq(0,1,.01)[34])] = "medium"
meta_data$Progenitor_sim[ prog_sim > quantile(prog_sim, seq(0,1,.01)[67])] = "high"
meta_data[not_sig_samples,"Progenitor_sim"] = "not_sig"

hsc_sim = log(res_coeff[,"hsc"]+1)
meta_data$HSC_sim = rep("low", length(hsc_sim))
meta_data$HSC_sim[ hsc_sim > quantile(hsc_sim, seq(0,1,.01)[34])] = "medium"
meta_data$HSC_sim[ hsc_sim > quantile(hsc_sim, seq(0,1,.01)[67])] = "high"
meta_data[not_sig_samples,"HSC_sim"] = "not_sig"

###

#absolute_mat = cbind(rowSums(res_coeff[,c("alpha","beta","gamma","delta")]), res_coeff[,"e13.5"],res_coeff[,"hsc"])

maxi = apply( res_coeff , FUN = which.max, MARGIN = 1 )
meta_data$Diff_Type = colnames(res_coeff)[maxi]
meta_data[not_sig_samples,"Diff_Type"] = "not_sig"
rownames(meta_data) = meta_data$Name
