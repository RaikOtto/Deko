
### RUN VARIANCE SELECTION FIRST

B_differentiated_model = readRDS(model_path)
B_differentiated = B_differentiated_model[[1]]
max_vec = B_differentiated_model[[2]]
#B_differentiated_six = readRDS("~/ArtDeco/inst/Models/Six_differentiation_stages.RDS")
B_differentiated_progenitor = readRDS("~/ArtDeco/inst/Models/Four_differentiation_stages_Lawlor_Progenitor.RDS")
B_undifferentiated = readRDS("~/ArtDeco/inst/Models/Four_differentiation_stages_Segerstolpe_Prog_Hisc.RDS")
#B_undifferentiated = readRDS("~/ArtDeco/inst/Models/Four_differentiation_stages_Baron_Prog_Stanescu_Hisc_Haber.RDS")

nr_permutations = 100
fit_differentiated = bseqsc_proportions(eset, B_differentiated, verbose = FALSE, absolute = T, log = F, perm = nr_permutations)
#fit_differentiated_six = bseqsc_proportions(eset, B_differentiated_six, verbose = FALSE, absolute = T, log = F, perm = nr_permutations)
fit_progenitor = bseqsc_proportions(eset, B_differentiated_progenitor, verbose = FALSE, absolute = T, log = F, perm = nr_permutations)
fit_undifferentiated = bseqsc_proportions(eset, B_undifferentiated, verbose = FALSE, absolute = T, log = F, perm = nr_permutations)

names(Marker_Gene_List) = str_to_lower(names(Marker_Gene_List))
eislet = new("ExpressionSet", exprs = as.matrix(expression_training_mat))

meta_data = meta_info[colnames(expression_training_mat),]
dim(meta_data)

sub_list = str_to_lower(meta_data$Subtype)
names(sub_list) = names(meta_data$Subtype)
fData(eislet) = data.frame( sub_list  )
pData(eislet) = data.frame( sub_list )
names(Marker_Gene_List)[!(names(Marker_Gene_List) %in% sub_list)]

B = bseqsc_basis(
    eislet,
    Marker_Gene_List,
    clusters = 'sub_list',
    samples = colnames(exprs(eislet)),
    ct.scale = FALSE
)
plotBasis(B, Marker_Gene_List, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')

### RUN VARIANCE SELECTION FIRST

eset = new("ExpressionSet", exprs=as.matrix(bam_data));

fit = bseqsc_proportions(eset, B, verbose = TRUE, absolute = T, log = F, perm = 100)

source("~/Deko/Scripts/Utility_script.R")
#meta_data$Diff_Type =  c("Differentiated","Progenitor","HSC")[maxi]
