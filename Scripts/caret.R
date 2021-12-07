library(mltest)
library(stringr)

<<<<<<< HEAD
<<<<<<< HEAD
performance_mat_deconvolution_grading_tertiary = read.table("~/Dropbox/testproject/Results/PanNEN_only/Deconvolution.Grading_tertiary.tsv",header = TRUE,  as.is = TRUE)
ml_predictions = ml_test(performance_mat_deconvolution_grading_tertiary[,2], performance_mat_deconvolution_grading_tertiary[,1], output.as.table = TRUE)
ml_predictions = round(ml_predictions,2)
write.table(ml_predictions,"~/Dropbox/testproject/Results/PanNEN_only/Visualization/Deconvolution.Grading_tertiary.tsv", sep ="\t", quote =F)

###

performance_mat_expression_grading_tertiary = read.table("~/Dropbox/testproject/Results/PanNEN_only/Expression_grading_tertiary_6_studies.tsv",header = TRUE,  as.is = TRUE)
ml_predictions = ml_test(performance_mat_expression_grading_tertiary[,2], performance_mat_expression_grading_tertiary[,1], output.as.table = TRUE)
ml_predictions = round(ml_predictions,2)
write.table(ml_predictions,"~/Dropbox/testproject/Results/PanNEN_only/Visualization/Expression_grading_tertiary_6_studies.tsv", sep ="\t", quote =F)

###

performance_mat_deconvolution_grading_binary = read.table("~/Dropbox/testproject/Results/PanNEN_only/Deconvolution.6_studies.hold_out.Grading_binary.tsv",header = TRUE,  as.is = TRUE)
ml_predictions = ml_test(performance_mat_deconvolution_grading_binary[,2], performance_mat_deconvolution_grading_binary[,1], output.as.table = TRUE)
ml_predictions = round(ml_predictions,2)
write.table(ml_predictions,"~/Dropbox/testproject/Results/PanNEN_only/Visualization/Deconvolution.6_studies.hold_out.Grading_binary.tsv", sep ="\t", quote =F)

#

performance_mat_expression_grading_binary = read.table("~/Dropbox/testproject/Results/PanNEN_only/Expression_grading_binary.3471_genes.hold_out.tsv",header = TRUE,  as.is = TRUE)
ml_predictions = ml_test(performance_mat_expression_grading_binary[,2], performance_mat_expression_grading_binary[,1], output.as.table = TRUE)
ml_predictions = round(ml_predictions,2)
write.table(ml_predictions,"~/Dropbox/testproject/Results/PanNEN_only/Visualization/Expression_grading_binary.3471_genes.hold_out.tsv", sep ="\t", quote =F)

###

performance_mat_deconvolution_net_nec = read.table("~/Dropbox/testproject/Results/PanNEN_only/Deconvolution.NET_NEC.6_studies.hold_out.tsv",header = TRUE,  as.is = TRUE)
ml_predictions = ml_test(performance_mat_deconvolution_net_nec[,2],true =  performance_mat_deconvolution_net_nec[,1], output.as.table = TRUE)
ml_predictions = round(ml_predictions,2)
write.table(ml_predictions,"~/Dropbox/testproject/Results/PanNEN_only/Visualization/Deconvolution.NET_NEC.6_studies.hold_out.tsv", sep ="\t", quote =F)

#

performance_mat_expression_net_nec = read.table("~/Dropbox/testproject/Results/PanNEN_only/Expression_NET_NEC_6_studies_3471_genes_hold_out.tsv",header = TRUE,  as.is = TRUE)
ml_predictions = ml_test(performance_mat_expression_net_nec[,2], performance_mat_expression_net_nec[,1], output.as.table = TRUE)
ml_predictions = round(ml_predictions,2)
write.table(ml_predictions,"~/Dropbox/testproject/Results/PanNEN_only/Visualization/Expression_NET_NEC_6_studies_3471_genes_hold_out.tsv", sep ="\t", quote =F)

=======
performance_mat_deconvolution_grading_tertiary = read.table("~/Downloads/Deconvolution.Grading_tertiary.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_deconvolution[,2], performance_mat_deconvolution[,1], output.as.table = TRUE)
classifier_metrics
=======
performance_mat_deconvolution_grading_tertiary = read.table("~/Dropbox/testproject/Results/PanNEN_only/Deconvolution.Grading_tertiary.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_deconvolution_grading_tertiary[,2], performance_mat_deconvolution_grading_tertiary[,1], output.as.table = TRUE)
write.table(classifier_metrics, "~/Dropbox/testproject/Results/PanNEN_only/Visualization/Deconvolution.Grading_tertiary.tsv")
>>>>>>> 031963ef0f9262657dafcbb781fe690c1f5ddc45

###

performance_mat_expression_grading_tertiary = read.table("~/Dropbox/testproject/Results/PanNEN_only/Expression_grading_tertiary_6_studies.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_expression_grading_tertiary[,2], performance_mat_expression_grading_tertiary[,1], output.as.table = TRUE)
write.table(classifier_metrics, "~/Dropbox/testproject/Results/PanNEN_only/Visualization/Expression_grading_tertiary_6_studies.tsv", sep ="\t", quote =F)

###

performance_mat_deconvolution_grading_binary = read.table("~/Dropbox/testproject/Results/PanNEN_only/Deconvolution.Grading_binary.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_deconvolution_grading_binary[,2], performance_mat_deconvolution_grading_binary[,1], output.as.table = TRUE)
write.table(classifier_metrics, "~/Dropbox/testproject/Results/PanNEN_only/Visualization/Deconvolution.6_studies.hold_out.Grading_binary.tsv", sep ="\t", quote =F)

###

performance_mat_expression_grading_binary = read.table("~/Dropbox/testproject/Results/PanNEN_only/Expression_grading_binary_6_studies.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_expression_grading_binary[,2], performance_mat_expression_grading_binary[,1], output.as.table = TRUE)
write.table(classifier_metrics, "~/Dropbox/testproject/Results/PanNEN_only/Visualization/Expression_grading_binary.3471_genes.hold_out.tsv", sep ="\t", quote =F)

###

performance_mat_deconvolution_NET_NEC = read.table("~/Dropbox/testproject/Results/PanNEN_only/Deconvolution.NET_NEC.6_studies.hold_out.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_deconvolution_NET_NEC[,2], performance_mat_deconvolution_NET_NEC[,1], output.as.table = TRUE)
write.table(classifier_metrics, "~/Dropbox/testproject/Results/PanNEN_only/Visualization/Deconvolution.NET_NEC.6_studies.hold_out.tsv", sep ="\t", quote =F)

###

performance_mat_expression_NET_NEC = read.table("~/Dropbox/testproject/Results/PanNEN_only/Expression_NET_NEC_6_studies_3471_genes_hold_out.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_expression_NET_NEC[,2], performance_mat_expression_NET_NEC[,1], output.as.table = TRUE)
<<<<<<< HEAD
classifier_metrics
>>>>>>> 202dfd25bed0797af48fb60fbfaa528be49ef295
=======
write.table(classifier_metrics, "~/Dropbox/testproject/Results/PanNEN_only/Visualization/Expression_NET_NEC_6_studies_3471_genes_hold_out.tsv", sep ="\t", quote =F)

###

single_performance_indices = read.table("~/Dropbox/testproject/Results/PanNEN_only/Visualization/ML_predictions_means.tsv",  sep ="\t", header = T)
aggregate(single_performance_indices[single_performance_indices$Type == "Deconvolution","Accuracy"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)
aggregate(single_performance_indices[single_performance_indices$Type == "Deconvolution","F1"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)
aggregate(single_performance_indices[single_performance_indices$Type == "Deconvolution","PPV"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)
aggregate(single_performance_indices[single_performance_indices$Type == "Deconvolution","Sensitivity"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)
aggregate(single_performance_indices[single_performance_indices$Type == "Deconvolution","Specificity"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)

aggregate(single_performance_indices[single_performance_indices$Type == "Baseline","Accuracy"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)
aggregate(single_performance_indices[single_performance_indices$Type == "Baseline","F1"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)
aggregate(single_performance_indices[single_performance_indices$Type == "Baseline","PPV"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)
aggregate(single_performance_indices[single_performance_indices$Type == "Baseline","Sensitivity"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)
aggregate(single_performance_indices[single_performance_indices$Type == "Baseline","Specificity"], by = list(single_performance_indices[single_performance_indices$Type == "Deconvolution","Experiment"]), FUN =mean)
>>>>>>> 031963ef0f9262657dafcbb781fe690c1f5ddc45
