aka3 = list(
    Histology   = c(
        Pancreatic_NEN = "BLACK",
        Colorectal_NEN = "Orange",
        Small_intestinal_NEN = "Yellow",
        Gastric_NEN = "purple",
        Liver = "Darkgreen",
        CUP = "pink"),
    Deco_type = c(alpha = "Blue",beta = "Yellow",gamma = "Orange",delta = "Purple", e13.5 = "Brown", hsc = "white"),
    NEC_NET = c(NEC= "red", NET = "blue", "NA" = "white"),
    Study = c(Groetzinger = "brown", Scarpa = "darkgreen"),
    MKI67 = c(high = "White", medium = "gray", low = "black"),
    INS = c(high = "White", low = "yellow"),
    Correlation = c(high = "Red", medium = "Yellow", low = "Green"),
    Deco_similarity = c(low = "Red", medium = "Yellow", high = "Green"),
    ImmuneScore = c(high = "White", medium = "gray", low = "Black"),
    StromalScore = c(high = "White", medium = "gray", low = "Black"),
    #Alpha_sim = c(low = "white", medium = "yellow", high = "blue",not_sig = "gray"),
    #Beta_sim = c(low = "white", medium = "yellow", high = "Yellow",not_sig = "gray"),
    #Gamma_sim = c(low = "white", medium = "yellow", high = "orange",not_sig = "gray"),
    #Delta_sim = c(low = "white", medium = "yellow", high = "Purple",not_sig = "gray"),
    #Ductal_sim = c(low = "white", medium = "yellow", high = "Black",not_sig = "gray"),
    #Acinar_sim = c(low = "white", medium = "yellow", high = "Brown",not_sig = "gray"),
    Alpha_sim = c(low = "white", medium = "yellow", high = "blue",not_sig = "gray"),
    Beta_sim = c(low = "white", medium = "Yellow", high = "darkgreen",not_sig = "gray"),
    Gamma_sim = c(low = "white", medium = "yellow", high = "orange",not_sig = "gray"),
    Delta_sim = c(low = "white", medium = "yellow", high = "Purple",not_sig = "gray"),
    Ductal_sim = c(low = "white", medium = "yellow", high = "Black",not_sig = "gray"),
    Acinar_sim = c(low = "white", medium = "yellow", high = "Brown",not_sig = "gray"),
    Progenitor_sim = c(low = "white", medium = "yellow", high = "red",not_sig = "gray"),
    Differentiated_sim_three= c(low = "white", medium = "yellow", high = "darkgreen",not_sig = "gray"),
    Diff_Type_Three = c(Differentiated = "darkgreen", Progenitor = "Brown", HESC = "Black"),
    Hsc_sim = c(low = "white", medium = "yellow", high = "black",not_sig = "gray"),
    Marker_Genes = c(high = "White", medium = "Yellow", low = "Red"),
    Functionality = c( Unknown = "White",Functional = "green", Non_Functional="red"),
    Grading = c( G1 = "Green",G2 = "Yellow", G3 = "Red", G0 = "white"),
    Subtype_Sadanandam = c(Norm = "darkgreen", Insulinoma = "Blue", MLP = "Orange", Intermediate = "brown", Unknown = "White"),
    #Differentiation_type = c( Differentiated = "green", Progenitor = "orange", HESC = "brown", not_sig = "gray")
    Differentiation_type = c( Alpha = "blue", Beta ="green", Gamma = "Orange", Delta = "Purple", Progenitor = "Yellow", HESC = "brown", not_sig = "gray")
)

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
    
        alpha_sim_scalar            = log( res_coeff[,"alpha"] + 1 )
        meta_data$Alpha_sim_scalar  = rep("",length(alpha_sim_scalar))
        meta_data$Alpha_sim         = rep("",length(alpha_sim_scalar))
        meta_data$Alpha_sim_scalar  = alpha_sim_scalar
        meta_data$Alpha_sim         = alpha_sim_scalar
        meta_data$Alpha_sim[ alpha_sim <= quantile(alpha_sim_scalar, seq(0,1,.01)[34])] = "low"
        meta_data$Alpha_sim[ alpha_sim > quantile(alpha_sim_scalar, seq(0,1,.01)[34])] = "medium"
        meta_data$Alpha_sim[ alpha_sim > quantile(alpha_sim_scalar, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Alpha_sim"] = "not_sig"
        
        Beta_sim_scalar            = log( res_coeff[,"beta"] + 1 )
        meta_data$Beta_sim_scalar  = rep("",length(Beta_sim_scalar))
        meta_data$Beta_sim         = rep("",length(Beta_sim_scalar))
        meta_data$Beta_sim_scalar  = Beta_sim_scalar
        meta_data$Beta_sim         = Beta_sim_scalar
        meta_data$Beta_sim[ Beta_sim_scalar <= quantile(Beta_sim_scalar, seq(0,1,.01)[34])] = "low"
        meta_data$Beta_sim[ Beta_sim_scalar > quantile(Beta_sim_scalar, seq(0,1,.01)[34])] = "medium"
        meta_data$Beta_sim[ Beta_sim_scalar > quantile(Beta_sim_scalar, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Beta_sim"] = "not_sig"
        
        Gamma_sim_scalar            = log( res_coeff[,"gamma"] + 1 )
        meta_data$Gamma_sim_scalar  = rep("",length(Gamma_sim_scalar))
        meta_data$Gamma_sim         = rep("",length(Gamma_sim_scalar))
        meta_data$Gamma_sim_scalar  = Gamma_sim_scalar
        meta_data$Gamma_sim         = Gamma_sim_scalar
        meta_data$Gamma_sim[ Gamma_sim_scalar <= quantile(Gamma_sim_scalar, seq(0,1,.01)[34])] = "low"
        meta_data$Gamma_sim[ Gamma_sim_scalar > quantile(Gamma_sim_scalar, seq(0,1,.01)[34])] = "medium"
        meta_data$Gamma_sim[ Gamma_sim_scalar > quantile(Gamma_sim_scalar, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Gamma_sim"] = "not_sig"
        
        Delta_sim_scalar            = log( res_coeff[,"delta"] + 1 )
        meta_data$Delta_sim_scalar  = rep("",length(Delta_sim_scalar))
        meta_data$Delta_sim         = rep("",length(Delta_sim_scalar))
        meta_data$Delta_sim_scalar  = Delta_sim_scalar
        meta_data$Delta_sim         = Delta_sim_scalar
        meta_data$Delta_sim[ Delta_sim_scalar <= quantile(Delta_sim_scalar, seq(0,1,.01)[34])] = "low"
        meta_data$Delta_sim[ Delta_sim_scalar > quantile(Delta_sim_scalar, seq(0,1,.01)[34])] = "medium"
        meta_data$Delta_sim[ Delta_sim_scalar > quantile(Delta_sim_scalar, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Delta_sim"] = "not_sig"
        
        Differentiated_sim_scalar = alpha_sim_scalar + beta_sim_scalar + gamma_sim_scalar + delta_sim_scalar
        meta_data$Differentiated_sim_scalar     = rep("",length(Differentiated_sim_scalar))
        meta_data$Differentiated_sim            = rep("",length(Differentiated_sim_scalar))
        meta_data$Differentiated_sim_scalar     = Differentiated_sim_scalar
        meta_data$Differentiated_sim            = Differentiated_sim_scalar
        meta_data$Differentiated_sim[ Differentiated_sim_scalar <= quantile(Differentiated_sim_scalar, seq(0,1,.01)[34])] = "low"
        meta_data$Differentiated_sim[ Differentiated_sim_scalar > quantile(Differentiated_sim_scalar, seq(0,1,.01)[34])] = "medium"
        meta_data$Differentiated_sim[ Differentiated_sim_scalar > quantile(Differentiated_sim_scalar, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Differentiated_sim"] = "not_sig"
    
    }
    
    if ( nr_fit == 2 ){
    
        Acinar_sim_scalar = log( res_coeff[,"acinar"] + 1 )
        meta_data$Acinar_sim_scalar     = rep("",length(Acinar_sim_scalar))
        meta_data$Acinar_sim            = rep("",length(Acinar_sim_scalar))
        meta_data$Acinar_sim_scalar     = Acinar_sim_scalar
        meta_data$Acinar_sim            = Acinar_sim_scalar
        meta_data$Acinar_sim[ Acinar_sim_scalar <= quantile(Acinar_sim_scalar, seq(0,1,.01)[34])] = "low"
        meta_data$Acinar_sim[ Acinar_sim_scalar > quantile(Acinar_sim_scalar, seq(0,1,.01)[34])] = "medium"
        meta_data$Acinar_sim[ Acinar_sim_scalar > quantile(Acinar_sim_scalar, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Acinar_sim"] = "not_sig"
        
        Ductal_sim_scalar = log( res_coeff[,"ductal"] + 1 )
        meta_data$Ductal_sim_scalar     = rep("",length(Ductal_sim_scalar))
        meta_data$Ductal_sim            = rep("",length(Ductal_sim_scalar))
        meta_data$Ductal_sim_scalar     = Ductal_sim_scalar
        meta_data$Ductal_sim            = Ductal_sim_scalar
        meta_data$Ductal_sim[ Ductal_sim_scalar <= quantile(Ductal_sim_scalar, seq(0,1,.01)[34])] = "low"
        meta_data$Ductal_sim[ Ductal_sim_scalar > quantile(Ductal_sim_scalar, seq(0,1,.01)[34])] = "medium"
        meta_data$Ductal_sim[ Ductal_sim_scalar > quantile(Ductal_sim_scalar, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Ductal_sim"] = "not_sig"
    }
    
    if ( nr_fit == 3 ){
    
        Prog_sim_scalar            = log( res_coeff[,"pancreatic_progenitor"] + 1 )
        meta_data$Prog_sim_scalar  = rep("",length(Prog_sim_scalar))
        meta_data$Prog_sim         = rep("",length(Prog_sim_scalar))
        meta_data$Prog_sim_scalar  = Prog_sim_scalar
        meta_data$Prog_sim         = Prog_sim_scalar
        meta_data$Prog_sim[ Prog_sim_scalar <= quantile(Prog_sim_scalar, seq(0,1,.01)[34])] = "low"
        meta_data$Prog_sim[ Prog_sim_scalar > quantile(Prog_sim_scalar, seq(0,1,.01)[34])] = "medium"
        meta_data$Prog_sim[ Prog_sim_scalar > quantile(Prog_sim_scalar, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Prog_sim"] = "not_sig"
        
        Hsc_sim_scalar            = log( res_coeff[,"hsc"] + 1 )
        meta_data$Hsc_sim_scalar  = rep("",length(Hsc_sim_scalar))
        meta_data$Hsc_sim         = rep("",length(Hsc_sim_scalar))
        meta_data$Hsc_sim_scalar  = Hsc_sim_scalar
        meta_data$Hsc_sim         = Hsc_sim_scalar
        meta_data$Hsc_sim[ Hsc_sim_scalar <= quantile(Hsc_sim_scalar, seq(0,1,.01)[34])] = "low"
        meta_data$Hsc_sim[ Hsc_sim_scalar > quantile(Hsc_sim_scalar, seq(0,1,.01)[34])] = "medium"
        meta_data$Hsc_sim[ Hsc_sim_scalar > quantile(Hsc_sim_scalar, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Hsc_sim"] = "not_sig"

    }
}

res_coeff = meta_data[,c("Alpha_sim_scalar","Beta_sim_scalar","Gamma_sim_scalar","Delta_sim_scalar","Prog_sim_scalar","Hsc_sim_scalar")]
colnames(res_coeff) = c("Alpha","Beta", "Gamma","Delta","Progenitor","HESC")
p_value_mat = meta_data[,c("Alpha_sim","Beta_sim","Gamma_sim","Delta_sim","Prog_sim","Hsc_sim")]

maxi = apply( res_coeff , FUN = which.max, MARGIN = 1 )
meta_data$Differentiation_type = colnames(res_coeff)[maxi]
not_sig_samples = p_value_mat[meta_data$Name,]
for(i in 1:length(meta_data$Differentiation_type)){
    if( p_value_mat[i,maxi[i]] == "not_sig" )
        meta_data[i,"Differentiation_type"] = "not_sig"
}
