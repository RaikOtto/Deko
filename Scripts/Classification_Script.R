### RUN VARIANCE SELECTION FIRST

B_differentiated_model = readRDS("~/ArtDeco/inst/Models/Four_differentiation_stages_Lawlor.RDS")
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

fits = list(fit_differentiated, fit_progenitor, fit_undifferentiated)
meta_data <<- meta_info[colnames(eset),]

source("~/Deko/Scripts/Utility_script.R")
