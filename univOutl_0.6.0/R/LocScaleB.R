LocScaleB <- function(x, k = 3, method = "MAD", weights = NULL, id = NULL,
                      exclude = NA, logt = FALSE, return.dataframe = FALSE) {
    
    # ---- input checks --------------------------------------------------------
    method <- match.arg(tolower(method),
                        c("mad", "iqr", "dq", "idr", "dd",
                          "gini", "scaletau2", "qn", "sn"))
    
    n <- length(x)
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    if (!is.numeric(k) || length(k) != 1 || is.na(k) || k <= 0)
        stop("'k' must be a single positive number")
    
    if (is.null(id)) id <- seq_len(n)
    if (length(id) != n) stop("'id' must have the same length as 'x'")
    
    if (!is.null(weights)) {
        if (length(weights) != n)
            stop("'weights' must have the same length as 'x'")
        if (anyNA(weights) || any(weights < 0))
            stop("'weights' must be non-negative and cannot contain NAs")
        if (method %in% c("gini", "scaletau2", "qn", "sn"))
            stop("Weights cannot be used with method '", method, "'")
    }
    
    # ---- remove values in 'exclude' and NAs ----------------------------------
    tst  <- (x %in% exclude) | is.na(x)
    keep <- !tst
    
    yy  <- x[keep]
    lab <- id[keep]
    ww  <- weights[keep]          # NULL stays NULL
    
    if (length(yy) < 3)
        stop("At least 3 valid observations are needed")
    
    # ---- log transform -------------------------------------------------------
    if (logt) {
        if (any(yy <= -1)) stop("log(x+1) cannot be applied: values <= -1 found")
        yy <- log(yy + 1)
        warning("Please note that log(x+1) transformation is applied")
    }
    
    # ---- quartiles -----------------------------------------------------------
    if (is.null(weights)) {
        qq <- quantile(x = yy, probs = c(0.25, 0.50, 0.75))
    } else {
        qq <- Hmisc::wtd.quantile(x = yy, weights = ww, probs = c(0.25, 0.50, 0.75))
    }
    qq  <- unname(qq)
    med <- qq[2]
    
    if (all(abs(qq) < 1e-06)) stop("Quartiles are all equal to 0")
    
    # ---- scale measures ------------------------------------------------------
    if (method %in% c("idr", "dd")) {
        if (is.null(weights)) {
            qq10 <- quantile(x = yy, probs = c(0.10, 0.50, 0.90))
        } else {
            qq10 <- Hmisc::wtd.quantile(x = yy, weights = ww,
                                        probs = c(0.10, 0.50, 0.90))
        }
        qq10 <- unname(qq10)
        
        if (method == "idr") {
            ddl <- ddr <- (qq10[3] - qq10[1]) / 2.5631
        } else {
            bow <- ((qq10[3] - qq10[2]) - (qq10[2] - qq10[1])) / (qq10[3] - qq10[1])
            message("Bowley's type skewness coefficient is: ", round(bow, 4))
            ddl <- (qq10[2] - qq10[1]) / 1.2816
            ddr <- (qq10[3] - qq10[2]) / 1.2816
        }
    }
    
    if (method == "iqr") {
        ddl <- ddr <- (qq[3] - qq[1]) / 1.3490
    }
    
    if (method == "dq") {
        bow <- ((qq[3] - qq[2]) - (qq[2] - qq[1])) / (qq[3] - qq[1])
        message("Bowley's skewness coefficient is: ", round(bow, 4))
        ddl <- (qq[2] - qq[1]) / 0.6745
        ddr <- (qq[3] - qq[2]) / 0.6745
    }
    
    if (method == "mad") {
        if (is.null(weights)) {
            ddl <- ddr <- mad(yy)
        } else {
            absdev <- abs(yy - med)
            ddl <- ddr <- unname(Hmisc::wtd.quantile(x = absdev, weights = ww,
                                                     probs = 0.5)) * 1.4826
        }
    }
    
    if (method == "gini")      ddl <- ddr <- Hmisc::GiniMd(yy) / 2 * sqrt(pi)
    if (method == "scaletau2") ddl <- ddr <- robustbase::scaleTau2(yy)
    if (method == "qn")        ddl <- ddr <- robustbase::Qn(yy)
    if (method == "sn")        ddl <- ddr <- robustbase::Sn(yy)
    
    ddl <- unname(ddl)
    ddr <- unname(ddr)
    
    # ---- check scale estimates -----------------------------------------------
    if (is.na(ddl) || ddl <= 0) {
        warning("Estimated left scale measure is 0 or NA: ",
                "it will be substituted with |0.05*Median|")
        ddl <- abs(0.05 * med)
    }
    if (is.na(ddr) || ddr <= 0) {
        warning("Estimated right scale measure is 0 or NA: ",
                "it will be substituted with |0.05*Median|")
        ddr <- abs(0.05 * med)
    }
    if (ddl <= 0 || ddr <= 0)
        stop("Scale estimate is 0 even after the |0.05*Median| substitution")
    
    # ---- scores --------------------------------------------------------------
    ee <- yy - med
    zz <- ifelse(yy < med, ee / ddl, ee / ddr)
    
    # ---- bounds --------------------------------------------------------------
    low.b <- med - k * ddl
    up.b  <- med + k * ddr
    bounds <- c(lower = low.b, upper = up.b)
    
    # ---- outliers (by tail) ----------------------------------------------------
    lower <- (yy < low.b)
    upper <- (yy > up.b)
    outl  <- lower | upper
    
    if (!any(outl)) {
        message("No outliers found")
    } else {
        message("No. of outliers in left tail: ", sum(lower))
        message("No. of outliers in right tail: ", sum(upper), "\n")
    }
    
    # ---- output ----------------------------------------------------------------
    if (ddl == ddr) {
        pp <- c(scale = ddr)
    } else {
        pp <- c(sc.left = ddl, sc.right = ddr)
    }
    pp <- c(median = med, pp)
    
    fine0 <- list(pars = pp, bounds = bounds)
    
    if (return.dataframe) {
        to.check <- data.frame(id = id[tst], x = x[tst])
        if (!is.null(weights)) to.check$weight <- weights[tst]
        
        df.out <- data.frame(id = lab, x = x[keep])
        if (logt) df.out$log.x <- yy
        if (!is.null(weights)) df.out$weight <- ww
        df.out$score    <- zz
        df.out$outliers <- as.integer(outl)
        
        fine1 <- list(excluded = to.check, data = df.out)
    } else {
        fine1 <- list(excluded = id[tst],
                      outliers = lab[outl],
                      lowOutl  = lab[lower],
                      upOutl   = lab[upper])
    }
    
    c(fine0, fine1)
}