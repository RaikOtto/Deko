library("stringr")
library("limma")

transcriptome_files
i = 7
transcriptome_files[i]

dataset_query = tail(as.character(unlist(str_split(transcriptome_files[i],pattern = "/"))),1)
dataset_query = str_replace_all(dataset_query,".tsv","")

if (type == "ductal") {
    models = models_ductal
} else if  (type == "hisc") {
    models = models_hisc
}

dataset_training = as.character(unlist(models[((i-1) %% 3) + 1]))

transcriptome_file = transcriptome_files[i]
visualization_file = visualization_files[i]

print(i)

if (file.exists(path_benchmark_files)){
    benchmark_results_t = read.table(path_benchmark_files,sep="\t",stringsAsFactors = F,header = T)
    
    if (nrow(benchmark_results_t) >= i)
        next(paste0("Skipping ",i))
}




algorithm = str_to_lower(algorithm)
transcriptome_data = read.table(transcriptome_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
colnames(transcriptome_data) = str_replace_all(colnames(transcriptome_data),pattern="^X","")
meta_data = meta_info[colnames(transcriptome_data),]
if( !nrow(meta_data) == ncol(transcriptome_data))
    stop("Inconclusive meta data and transcriptome file dimensions")
row_names = as.character(rownames(transcriptome_data))
col_names = as.character(colnames(transcriptome_data))
transcriptome_data = matrix(as.double(as.character(unlist(transcriptome_data))),ncol = length(col_names))
rownames(transcriptome_data) = row_names
colnames(transcriptome_data) = col_names

visualization_data = read.table(visualization_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
colnames(visualization_data) = str_replace_all(colnames(visualization_data),pattern="^X","")

name_query_data = tail(as.character(unlist(str_split(transcriptome_file,pattern ="/"))),1)
name_query_data = str_replace_all(name_query_data,pattern =".tsv","")

name_training_data = tail(as.character(unlist(str_split(dataset_training,pattern ="_"))),1)

name_algorithm_datatype = tail(as.character(unlist(str_split(path_benchmark_files,pattern ="/"))),1)
name_algorithm_datatype = str_replace_all(name_algorithm_datatype,pattern ="(.tsv)|(Benchmark_results.)","")
name_algorithm = head(as.character(unlist(str_split(name_algorithm_datatype,pattern ="\\."))),1)
name_datatype = tail(as.character(unlist(str_split(name_algorithm_datatype,pattern ="\\."))),1)

identifier = tail(as.character(unlist(str_split(dataset_training[1],pattern = "_"))),1)
if (str_detect(dataset_training[2], "HISC")){
    dataset_training_label = paste(identifier,"hisc",sep="_")
} else {dataset_training_label = identifier}

ki_index = which(rownames(transcriptome_data) == "MKI67")

###

model_path = paste0(
    c(
        "~/Deko/Data/Bench_data/Models/",
        dataset_query,
        dataset_training[2],
        algorithm,
        "RDS"
    ),
    collapse = "."
)
model_path = str_replace_all(model_path, pattern = "/\\.","/")

###

if (! file.exists(model_path)){
    
    deconvolution_results = Determine_differentiation_stage(
        transcriptome_data = transcriptome_data,
        deconvolution_algorithm = str_to_lower(algorithm),
        models = dataset_training,
        nr_permutations = 1000,
        output_file = ""
    )
    saveRDS(deconvolution_results, file = model_path)
} else {
    deconvolution_results = readRDS(model_path)
}

if (sum( meta_data[rownames(deconvolution_results),"Grading"] != "") > 0){
    
    meta_data = meta_data[meta_data$Grading!="",]
    deconvolution_results = deconvolution_results[rownames(meta_data),]
    deconvolution_results$Grading = meta_data$Grading
    visualization_data  = visualization_data[,rownames(meta_data)]
    transcriptome_data  = transcriptome_data[,rownames(meta_data)]
    
}

if( length(ki_index) != 0 ){
    
    deconvolution_results[,"MKI67"] = rep(0,nrow(deconvolution_results))
    deconvolution_results[,"MKI67"] = log(as.double(transcriptome_data[ki_index[1],])+1)
    
} else {
    
    deconvolution_results[,"MKI67"] = as.double(meta_data$KI67)
}

if (str_detect( transcriptome_file,pattern = "Wiedenmann_Scarpa_GSE73338")){
    mki_67_wiedenmann = read.table("~/Deko/Data/Cancer_Pancreas_Bulk_Array/Wiedenmann_Scarpa/Groetzinger_Scarpa_57.primary_only.tsv",sep="\t",header = T)
    colnames(mki_67_wiedenmann) = str_replace_all(colnames(mki_67_wiedenmann),pattern="^X","")
    deconvolution_results[colnames(mki_67_wiedenmann),"MKI67"] = as.double(mki_67_wiedenmann["MKI67",])
}

deconvolution_results$Strength_de_differentiation[which(is.infinite(as.double(deconvolution_results$Strength_de_differentiation)))] = -4
deconvolution_results$Confidence_score_dif[which(is.infinite(as.double(deconvolution_results$Confidence_score_dif)))] = 0

vis_mat = create_visualization_matrix(
    visualization_data = visualization_data,
    deconvolution_results = deconvolution_results,
    confidence_threshold = confidence_threshold,
    high_threshold = high_threshold,
    low_threshold = low_threshold
)

### results parsing

graphics_path_heatmap = paste("~/Deko/Results/Images",algorithm, sep ="/")
graphics_path_heatmap = paste(graphics_path_heatmap,"Heatmap", sep ="/")
if (! dir.exists(graphics_path_heatmap))
    dir.create(graphics_path_heatmap)
graphics_path_heatmap = paste(graphics_path_heatmap,name_datatype, sep ="/")
if (! dir.exists(graphics_path_heatmap))
    dir.create(graphics_path_heatmap)
graphics_path_heatmap = paste(graphics_path_heatmap,name_query_data, sep ="/")
if (! dir.exists(graphics_path_heatmap))
    dir.create(graphics_path_heatmap)
graphics_path_heatmap = paste(graphics_path_heatmap,paste0(name_training_data,".pdf"),sep = "/")

pdf(graphics_path_heatmap,onefile = FALSE)#,width="1024px",height="768px")
create_heatmap_differentiation_stages(
    visualization_data = visualization_data,
    deconvolution_results = deconvolution_results,
    vis_mat = vis_mat,
    confidence_threshold = confidence_threshold,
    show_colnames = F,
    aggregate_differentiated_stages = FALSE,
    high_threshold = high_threshold,
    low_threshold = low_threshold
)
dev.off()

### survival curve

if (sum(meta_data$OS_Tissue != "") != 0){
    
    meta_data = meta_info[rownames(vis_mat),]
    vis_mat = vis_mat[rownames(deconvolution_results),]
    vis_mat$OS_Tissue = as.double(str_replace_all(meta_data$OS_Tissue, pattern = ",", "\\."))
    vis_mat$OS_Tissue[is.na(vis_mat$OS_Tissue)] = 1
    vis_mat$Grading = meta_data$Grading
    vis_mat$Zensur = meta_data$Zensur
    
    mki67 = deconvolution_results[rownames(vis_mat),"MKI67"]
    mki67[mki67 <= mean(mki67)] = "low"
    mki67[mki67 != "low"] = "high"
    ductal = deconvolution_results$ductal
    ductal[ductal <= mean(ductal)] = "low"
    ductal[ductal != "low"] = "high"
    ratio = vis_mat[rownames(deconvolution_results),"Ratio_numeric"]
    ratio[ratio <= mean(ratio)] = "low"
    ratio[ratio != "low"] = "high"
    
    #mki67 = vis_mat$MKI67
    #ductal = vis_mat$ductal
    #ratio = vis_mat$Ratio
    
    data = vis_mat[,c("OS_Tissue","Zensur")]
    
    graphics_path_survival = paste("~/Deko/Results/Images",algorithm, sep ="/")
    graphics_path_survival = paste(graphics_path_survival,"Survival", sep ="/")
    if (! dir.exists(graphics_path_survival))
        dir.create(graphics_path_survival)
    graphics_path_survival = paste(graphics_path_survival,name_datatype, sep ="/")
    if (! dir.exists(graphics_path_survival))
        dir.create(graphics_path_survival)
    graphics_path_survival = paste(graphics_path_survival,name_query_data, sep ="/")
    if (! dir.exists(graphics_path_survival))
        dir.create(graphics_path_survival)
    
    if ("hisc" %in% colnames(deconvolution_results)){
        
        hisc = deconvolution_results$hisc
        hisc[hisc <= mean(hisc)] = "low"
        hisc[hisc != "low"] = "high"
        
        #hisc = vis_mat$hisc
        
        data = cbind(mki67,ductal,hisc,ratio,data)
        
        surv_hisc = survminer::surv_pvalue(survival::survfit( survival::Surv( as.double(data$OS_Tissue) ) ~ hisc), data = data[data$hisc!="not_significant",], method = "survdiff")$pval
        if ( sum(hisc == "high") < 5) surv_hisc = 1
        
        graphics_path_survival_hisc = paste(graphics_path_survival,paste0(paste(name_training_data,"hisc",sep="_"),".pdf"),sep = "/")
        
        data = cbind(hisc,vis_mat[,c("OS_Tissue","Zensur")])
        data = data[data$hisc!="not_significant",]
        fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Zensur ) ~ data$hisc)
        
        pdf(graphics_path_survival_hisc,onefile = FALSE)#,width="1024px",height="768px")
        print(survminer::ggsurvplot(fit, data = data, risk.table = T, pval = T, censor.size = 10))
        dev.off()
        
    } else {
        
        data = cbind(mki67,ductal,ratio,data)
        surv_hisc = 1
    }
    
    data = data[ !is.na(data$OS_Tissue),]
    
    data$mki67 = mki67
    data$ductal = ductal
    data$ratio = ratio
    
    surv_mki67 = survminer::surv_pvalue(survival::survfit( survival::Surv( as.double(data$OS_Tissue) ) ~ data$mki67), data = data, method = "survdiff")$pval
    if ( sum(mki67 == "high") < 5) surv_mki67 = 1
    surv_ductal = survminer::surv_pvalue(survival::survfit( survival::Surv( as.double(data$OS_Tissue) ) ~ data$ductal), data = data, method = "survdiff")$pval
    if ( sum(ductal == "high") < 5) surv_ductal = 1
    surv_ratio = survminer::surv_pvalue(survival::survfit( survival::Surv( as.double(data$OS_Tissue) ) ~ data$ratio), data = data, method = "survdiff")$pval
    if ( sum(ratio == "high") < 5) surv_ratio = 1
    
    data = cbind(mki67,ductal,ratio,vis_mat[,c("OS_Tissue","Zensur")])
    
    # mki67
    
    graphics_path_survival_mki67 = paste(graphics_path_survival,paste0(paste(name_training_data,"mki67",sep="_"),".pdf"),sep = "/")
    
    fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Zensur ) ~ data$mki67)
    
    pdf(graphics_path_survival_mki67,onefile = FALSE)#,width="1024px",height="768px")
    print(survminer::ggsurvplot(fit, data = data, risk.table = F, pval = T, censor.size = 10))
    dev.off()
    
    # ductal 
    
    graphics_path_survival_ductal = paste(graphics_path_survival,paste0(paste(name_training_data,"ductal",sep="_"),".pdf"),sep = "/")
    
    data = cbind(ductal,vis_mat[,c("OS_Tissue","Zensur")])
    data = data[data$ductal != "not_significant",]
    fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Zensur ) ~ data$ductal)
    
    pdf(graphics_path_survival_ductal,onefile = FALSE)#,width="1024px",height="768px")
    print(survminer::ggsurvplot(fit, data = data, risk.table = F, pval = T, censor.size = 10))
    dev.off()
    
    # ratio
    
    graphics_path_survival_ratio = paste(graphics_path_survival,paste0(paste(name_training_data,"ratio",sep="_"),".pdf"),sep = "/")
    
    data = cbind(ratio,vis_mat[,c("OS_Tissue","Zensur")])
    fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Zensur ) ~ data$ratio)
    
    pdf(graphics_path_survival_ratio,onefile = FALSE)#,width="1024px",height="768px")
    print(survminer::ggsurvplot(fit, data = data, risk.table = F, pval = T, censor.size = 10))
    dev.off()
    
} else {
    surv_mki67 = surv_ductal = surv_hisc = surv_ratio = 1
}


###


## comparison deconvolution

# ductal

deco_grading = meta_info[rownames(deconvolution_results),"Grading"] %>%  str_replace_all(pattern = "G","") %>% as.double

cor.test(deco_grading, deconvolution_results$ductal)
cor.test(deco_grading, deconvolution_results$hisc)
cor.test(deco_grading, deconvolution_results$MKI67)
