library(bseqsc)

expression_training_mat = read.table(
    "~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Lawlor.Alpha.Beta.Gamma.Delta.tsv",
    sep ="\t",
    header = T,
    stringsAsFactors = F
)
colnames(expression_training_mat) = 
    str_replace(colnames(expression_training_mat), pattern = "\\.", "_")
colnames(expression_training_mat) =
    str_replace(colnames(expression_training_mat), pattern = "^X", "")

meta_data <<- meta_info[colnames(expression_training_mat),]

if (length(subtype_vector) == 0)
    stop(paste0("You have to provide the sample subtypes labels for model training"))
subtype_vector = str_to_lower(meta_data$Subtype)
Marker_Gene_List = list()

### Data cleansing

row_var = apply(expression_training_mat, FUN = var, MARGIN = 1)
col_var = apply(expression_training_mat, FUN = var, MARGIN = 2)
expression_training_mat = expression_training_mat[row_var != 0,col_var != 0]
expression_training_mat = expression_training_mat[rowSums(expression_training_mat) >= 1,]

training_nr_marker_genes = 100

source("~/Deko/Scripts/Marker_gene_identification.R")
for( subtype in unique(subtype_vector) ){
    Marker_Gene_List[[subtype]] = identify_marker_genes(
        expression_training_mat,
        subtype_vector,
        subtype,
        training_nr_marker_genes
    )
}
print("Finished extracting marker genes for subtypes")

# Prepare bseq training

training_mat_bseq = new(
    "ExpressionSet",
    exprs = as.matrix(bam_data_1)
)
fData(training_mat_bseq) = data.frame( subtype_vector )
pData(training_mat_bseq) = data.frame( subtype_vector )

Basis = bseqsc_basis(
    training_mat_bseq,
    Marker_Gene_List,
    clusters = 'subtype_vector',
    samples = colnames(training_mat_bseq),
    log = F,
    ct.scale = FALSE
)

print(
    "Basis trained, estimating deconvolution thresholds, this may take some time")

test_mat = new(
    "ExpressionSet",
    exprs = as.matrix(expression_training_mat)
);

training_nr_permutations = 100
fit = bseqsc_proportions(
    test_mat,
    Basis,
    verbose = FALSE,
    absolute = T,
    log = F,
    perm = training_nr_permutations
)

print("Finished threshold determination")

res_coeff = t(fit$coefficients)
res_coeff_mat = as.double(unlist(res_coeff))
res_coeff_mat = as.data.frame(
    matrix(
        res_coeff_mat,
        ncol = ncol(res_coeff),
        nrow = nrow(res_coeff)
    )
)
rownames(res_coeff_mat) = rownames(res_coeff)
colnames(res_coeff_mat) = colnames(res_coeff)
res_cor   = fit$stats

res_coeff[ is.na(res_coeff) ] = 0.0
res_cor[ is.na(res_cor) ] = 0.0

self_scores = list()
for (subtype in unique(subtype_vector)){
    self_scores[[subtype]] = as.double(
        res_coeff[
            which(subtype_vector == subtype),
            subtype
        ]
    )
}

model = list(Basis, self_scores, Marker_Gene_List)
model_path = "~/Deko/Models/Alpha_Beta_Gamma_Delta_Lawlor.RDS"
saveRDS(model,model_path)
