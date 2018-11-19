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
    Alpha_sim = c(low = "white", medium = "yellow", high = "blue"),
    Beta_sim = c(low = "white", medium = "Yellow", high = "darkgreen"),
    Gamma_sim = c(low = "white", medium = "yellow", high = "orange"),
    Delta_sim = c(low = "white", medium = "yellow", high = "Purple"),
    Ductal_sim = c(low = "white", medium = "yellow", high = "Black"),
    Acinar_sim = c(low = "white", medium = "yellow", high = "Brown"),
    Progenitor_sim = c(low = "white", medium = "yellow", high = "orange",not_sig = "gray"),
    Differentiated_sim = c(low = "white", medium = "yellow", high = "darkgreen",not_sig = "gray"),
    Differentiated_sim_three= c(low = "white", medium = "yellow", high = "darkgreen",not_sig = "gray"),
    Diff_Type_Three = c(Differentiated = "darkgreen", Progenitor = "Brown", HESC = "Black"),
    HESC_sim = c(low = "white", medium = "yellow", high = "darkred",not_sig = "gray"),
    Marker_Genes = c(high = "White", medium = "Yellow", low = "Red"),
    Functionality = c( Unknown = "White",Functional = "green", Non_Functional="red"),
    Grading = c( G1 = "Green",G2 = "Yellow", G3 = "Red", G0 = "white"),
    Subtype_Sadanandam = c(Norm = "darkgreen", Insulinoma = "Blue", MLP = "Orange", Intermediate = "brown", Unknown = "White"),
    Diff_Type = c( alpha = "blue", beta = "Green", gamma = "Orange", delta = "Purple"),
    Y_Subtype = c( Alpha = "blue", Beta = "Green", Gamma = "Orange", Delta = "Purple")
)

res_coeff = t(fit$coefficients)
res_cor   = fit$stats

res_coeff[ is.na(res_coeff) ] = 0.0
res_cor[ is.na(res_cor) ] = 0.0

meta_data = meta_info[ rownames(res_coeff),]

not_sig_samples = rownames(res_cor)[res_cor[,"P-value"] > .1]
not_sig_samples

#sig_quantiles = 

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

if (str_to_upper("ductal") %in% str_to_upper(names(pancreasMarkers))){

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

}

diff_sim = log(rowSums(res_coeff[,c("alpha","beta","gamma","delta")])+1)
meta_data$Differentiated_sim = rep("low", length(diff_sim))
meta_data$Differentiated_sim[ diff_sim > quantile(diff_sim, seq(0,1,.01)[34])] = "medium"
meta_data$Differentiated_sim[ diff_sim > quantile(diff_sim, seq(0,1,.01)[67])] = "high"
meta_data[not_sig_samples,"Differentiated_sim"] = "not_sig"

if (str_to_upper( "HESC") %in% str_to_upper(names(pancreasMarkers))){

        prog_sim = log(res_coeff[,"progenitor"]+1)
        meta_data$Progenitor_sim = rep("low", length(prog_sim))
        meta_data$Progenitor_sim[ prog_sim > quantile(prog_sim, seq(0,1,.01)[34])] = "medium"
        meta_data$Progenitor_sim[ prog_sim > quantile(prog_sim, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"Progenitor_sim"] = "not_sig"
        
        hesc_sim = log(res_coeff[,"hesc"]+1)
        meta_data$HESC_sim = rep("low", length(hesc_sim))
        meta_data$HESC_sim[ hesc_sim > quantile(hesc_sim, seq(0,1,.01)[34])] = "medium"
        meta_data$HESC_sim[ hesc_sim > quantile(hesc_sim, seq(0,1,.01)[67])] = "high"
        meta_data[not_sig_samples,"HESC_sim"] = "not_sig"

}
###

#absolute_mat = cbind(rowSums(res_coeff[,c("alpha","beta","gamma","delta")]), res_coeff[,"e13.5"],res_coeff[,"hsc"])

res_coeff_aggregated = cbind(rowSums(res_coeff[,1:4]), res_coeff[,5:ncol(res_coeff)])
colnames(res_coeff_aggregated) = c("Differentiated","Progenitor","HESC")

maxi = apply( res_coeff , FUN = which.max, MARGIN = 1 )
meta_data$Diff_Type = colnames(res_coeff)[maxi]
maxi_three = apply( res_coeff_aggregated , FUN = which.max, MARGIN = 1 )
meta_data$Diff_Type_Three = colnames(res_coeff_aggregated[meta_data$Name,])[maxi_three]
meta_data[not_sig_samples,"Diff_Type"] = "not_sig"
rownames(meta_data) = meta_data$Name

meta_data$Diff_Type[meta_data$Diff_Type == "alpha"] = "Alpha"
meta_data$Diff_Type[meta_data$Diff_Type == "beta"] = "Beta"
meta_data$Diff_Type[meta_data$Diff_Type == "gamma"] = "Gamma"
meta_data$Diff_Type[meta_data$Diff_Type == "delta"] = "Delta"
meta_data$Diff_Type[meta_data$Diff_Type == "ductal"] = "Ductal"
meta_data$Diff_Type[meta_data$Diff_Type == "acinar"] = "Acinar"

meta_data$Differentiated_sim_three = three_sim = log(res_coeff_aggregated[,"Differentiated"]+1)
meta_data$Differentiated_sim_three[ three_sim <= quantile(three_sim, seq(0,1,.01)[34])] = "low"
meta_data$Differentiated_sim_three[ three_sim > quantile(three_sim, seq(0,1,.01)[34])] = "medium"
meta_data$Differentiated_sim_three[ three_sim > quantile(three_sim, seq(0,1,.01)[67])] = "high"
meta_data[not_sig_samples,"Differentiated_sim_three"] = "not_sig"
