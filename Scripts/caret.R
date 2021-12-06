library(mltest)

performance_mat_deconvolution_grading_tertiary = read.table("~/Downloads/Deconvolution.Grading_tertiary.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_deconvolution[,2], performance_mat_deconvolution[,1], output.as.table = TRUE)
classifier_metrics

###

performance_mat_expression_grading_tertiary = read.table("~/Downloads/Expression_grading_tertiary.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_expression_grading_tertiary[,2], performance_mat_expression_grading_tertiary[,1], output.as.table = TRUE)
classifier_metrics

###

performance_mat_deconvolution_grading_binary = read.table("~/Downloads/Deconvolution.Grading_binary.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_deconvolution_grading_binary[,2], performance_mat_deconvolution_grading_binary[,1], output.as.table = TRUE)
classifier_metrics

###

performance_mat_expression_grading_binary = read.table("~/Downloads/Expression_grading_binary.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_expression_grading_binary[,2], performance_mat_expression_grading_binary[,1], output.as.table = TRUE)
classifier_metrics

###

performance_mat_deconvolution_NET_NEC = read.table("~/Downloads/Deconvolution.NET_NEC.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_deconvolution_NET_NEC[,2], performance_mat_deconvolution_NET_NEC[,1], output.as.table = TRUE)
classifier_metrics

###

performance_mat_expression_NET_NEC = read.table("~/Downloads/Expression_NET_NEC.tsv",header = TRUE,  as.is = TRUE)
classifier_metrics <- ml_test(performance_mat_expression_NET_NEC[,2], performance_mat_expression_NET_NEC[,1], output.as.table = TRUE)
classifier_metrics
