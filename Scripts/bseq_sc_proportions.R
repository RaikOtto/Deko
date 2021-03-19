bseqsc_proportions = function (x, reference = NULL, log = NULL, ..., verbose = TRUE) 
{
    message <- if (verbose) 
        message
    else function(...) invisible(NULL)
    if (is.null(reference)) {
        message("* Using pancreatic islet reference basis matrix: ", 
                appendLF = FALSE)
        reference <- ldata("PancreasIslet")
    }
    CIBERSORT <- bseqsc_config(error = TRUE)
    y <- x
    if (is(y, "ExpressionSet")) 
        y <- exprs(y)
    ids <- intersect(rownames(y), rownames(reference))
    message(sprintf("* Data features: %s", str_out(rownames(y), 
                                                   total = TRUE)))
    message(sprintf("* Basis features: %s", str_out(rownames(reference), 
                                                    total = TRUE)))
    message(sprintf("* Common features: %s", str_out(ids, total = TRUE)))
    islog <- log %||% is_logscale(y)
    if (islog) {
        message("* Converting to linear scale")
        y <- expb(y, 2)
    }
    tdir <- tempfile("CIBERSORT")
    dir.create(tdir)
    owd <- setwd(tdir)
    on.exit({
        unlink(tdir, recursive = TRUE)
        setwd(owd)
    })
    message("* Writing input files ... ", appendLF = FALSE)
    write.table(reference, file = xf <- "reference.tsv", sep = "\t", 
                row.names = TRUE, col.names = NA)
    write.table(y, file = yf <- "mixture.tsv", sep = "\t", row.names = TRUE, 
                col.names = NA)
    message("OK")
    message("* Running CIBERSORT ... ", appendLF = FALSE)
    res <- CIBERSORT(xf, yf, ...)
    message("OK")
    stats <- setdiff(colnames(res), colnames(reference))
    list(coefficients = t(res[, !colnames(res) %in% stats]), 
         stats = res[, stats])
}