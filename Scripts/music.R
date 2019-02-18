music_prop = function (
    bulk.eset,
    sc.eset,
    markers = NULL,
    clusters,
    samples,
    select.ct = NULL,
    ct.cov = FALSE,
    verbose = TRUE,
    iter.max = 1000,
    nu = 1e-04,
    eps = 0.01,
    inter.as.bseq = F,
    normalize = F
) 
{
    exprs(bulk.eset) = as.data.frame(exprs(bulk.eset))
    bulk.gene = rownames(bulk.eset)[rowMeans(exprs(bulk.eset)) != 0]
    bulk.eset = bulk.eset[bulk.gene, , drop = FALSE]
    if (is.null(markers)) {
        sc.markers = bulk.gene
    }
    else {
        sc.markers = intersect(bulk.gene, unlist(markers))
    }
    sc.basis = music_basis(sc.eset, non.zero = TRUE, markers = sc.markers, 
                           clusters = clusters, samples = samples, select.ct = select.ct, 
                           ct.cov = ct.cov, verbose = verbose)
    cm.gene = intersect(rownames(sc.basis$Disgn.mtx), bulk.gene)
    if (is.null(markers)) {
        if (length(cm.gene) < 0.2 * min(length(bulk.gene), nrow(sc.eset))) 
            stop("Too few common genes!")
    }
    else {
        if (length(cm.gene) < 0.2 * length(unlist(markers))) 
            stop("Too few common genes!")
    }
    if (verbose) {
        message(paste("Used", length(cm.gene), "common genes..."))
    }
    m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx))
    m.bulk = match(cm.gene, bulk.gene)
    D1 = sc.basis$Disgn.mtx[m.sc, ]
    M.S = colMeans(sc.basis$S, na.rm = T)
    Yjg = relative.ab(exprs(bulk.eset)[m.bulk, ])
    N.bulk = ncol(bulk.eset)
    if (ct.cov) {
        Sigma.ct = sc.basis$Sigma.ct[, m.sc]
        if (sum(Yjg[, i] == 0) > 0) {
            D1.temp = D1[Yjg[, i] != 0, ]
            Yjg.temp = Yjg[Yjg[, i] != 0, i]
            Sigma.ct.temp = Sigma.ct[, Yjg[, i] != 0]
            if (verbose) 
                message(paste(colnames(Yjg)[i], "has common genes", 
                              sum(Yjg[, i] != 0), "..."))
        }
        else {
            D1.temp = D1
            Yjg.temp = Yjg[, i]
            Sigma.ct.temp = Sigma.ct
            if (verbose) 
                message(paste(colnames(Yjg)[i], "has common genes", 
                              sum(Yjg[, i] != 0), "..."))
        }
        lm.D1.weighted = music.iter.ct(Yjg.temp, D1.temp, M.S, 
                                       Sigma.ct.temp, iter.max = iter.max, nu = nu, eps = eps, 
                                       inter.as.bseq = inter.as.bseq, normalize = normalize)
        Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
        Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
        weight.gene.temp = rep(NA, nrow(Yjg))
        weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
        Weight.gene = cbind(Weight.gene, weight.gene.temp)
        r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
        Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
    else {
        Sigma = sc.basis$Sigma[m.sc, ]
        valid.ct = (colSums(is.na(Sigma)) == 0) & (colSums(is.na(D1)) == 
                                                       0) & (!is.na(M.S))
        if (sum(valid.ct) <= 1) {
            stop("Not enough valid cell type!")
        }
        if (verbose) {
            message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
        }
        D1 = D1[, valid.ct]
        M.S = M.S[valid.ct]
        Sigma = Sigma[, valid.ct]
        Est.prop.allgene = NULL
        Est.prop.weighted = NULL
        Weight.gene = NULL
        r.squared.full = NULL
        Var.prop = NULL
        for (i in 1:N.bulk) {
            if (sum(Yjg[, i] == 0) > 0) {
                D1.temp = D1[Yjg[, i] != 0, ]
                Yjg.temp = Yjg[Yjg[, i] != 0, i]
                Sigma.temp = Sigma[Yjg[, i] != 0, ]
                if (verbose) 
                    message(paste(colnames(Yjg)[i], "has common genes", 
                                  sum(Yjg[, i] != 0), "..."))
            }
            else {
                D1.temp = D1
                Yjg.temp = Yjg[, i]
                Sigma.temp = Sigma
                if (verbose) 
                    message(paste(colnames(Yjg)[i], "has common genes", 
                                  sum(Yjg[, i] != 0), "..."))
            }
            lm.D1.weighted = music.iter(Yjg.temp, D1.temp, M.S, 
                                        Sigma.temp, iter.max = iter.max, nu = nu, eps = eps, 
                                        inter.as.bseq = inter.as.bseq, normalize = normalize)
            Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
            Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
            weight.gene.temp = rep(NA, nrow(Yjg))
            weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
            Weight.gene = cbind(Weight.gene, weight.gene.temp)
            r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
            Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
        }
    }
    colnames(Est.prop.weighted) = colnames(D1)
    rownames(Est.prop.weighted) = colnames(Yjg)
    colnames(Est.prop.allgene) = colnames(D1)
    rownames(Est.prop.allgene) = colnames(Yjg)
    names(r.squared.full) = colnames(Yjg)
    colnames(Weight.gene) = colnames(Yjg)
    rownames(Weight.gene) = cm.gene
    colnames(Var.prop) = colnames(D1)
    rownames(Var.prop) = colnames(Yjg)
    return(list(Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene, 
                Weight.gene = Weight.gene, r.squared.full = r.squared.full, 
                Var.prop = Var.prop))
}