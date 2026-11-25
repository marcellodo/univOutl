boxB <- function(x, k = 1.5,
                 method = c('asymmetric', 'asymmetric2', 'resistant',
                            'standard', 'adjbox'),
                 weights = NULL, id = NULL, exclude = NA, verbose = TRUE){

    # --- argument checks ---------------------------------------------------
    method <- match.arg(method)               # catches typos / invalid methods
    n <- length(x)
    if(n == 0) stop("'x' has length 0")
    if(!is.numeric(k) || length(k) != 1 || is.na(k) || k < 0)
        stop("'k' must be a single non-negative number")

    if(is.null(id)) id <- seq_len(n)          # safe for length 0
    if(length(id) != n) stop("'id' must have the same length as 'x'")

    if(!is.null(weights)){
        if(length(weights) != n) stop("'weights' must have the same length as 'x'")
        if(anyNA(weights) || any(weights < 0))
            stop("'weights' must be non-negative and without NA")
        # adjbox is fully unweighted (the medcouple has no weighted version)
        if(method == 'adjbox'){
            warning("method='adjbox' does not use weights: 'weights' ignored")
            weights <- NULL
        } else if(!requireNamespace("Hmisc", quietly = TRUE)){
            stop("Package 'Hmisc' is needed for weighted quantiles")
        }
    }
    if(method == 'adjbox' && !requireNamespace("robustbase", quietly = TRUE))
        stop("Package 'robustbase' is needed for method='adjbox'")

    # --- exclusions --------------------------------------------------------
    # NA/NaN are always excluded (they would make quantile() fail),
    # plus any value listed in 'exclude'
    tst <- (x %in% exclude) | is.na(x)
    to.check <- id[tst]
    yy  <- x[!tst]
    lab <- id[!tst]
    ww  <- if(is.null(weights)) NULL else weights[!tst]

    if(length(yy) == 0) stop("no values left after exclusions")

    # --- quartiles ---------------------------------------------------------
    if(is.null(weights)){
        qq <- quantile(x = yy, probs = c(0.25, 0.50, 0.75))
    } else {
        # note: weighted (type 'quantile') and unweighted (type 7) quartile
        # definitions differ slightly; normwt=FALSE -> weights as frequencies
        qq <- Hmisc::wtd.quantile(x = yy, weights = ww,
                                  probs = c(0.25, 0.50, 0.75))
    }

    # --- fences ------------------------------------------------------------
    if(method %in% c('resistant', 'standard')){
        vql <- vqu <- qq[3] - qq[1]
        if(abs(vql) < 1e-06) warning("IQR is 0")
        low.b <- qq[1] - k * vql
        up.b  <- qq[3] + k * vqu
    }

    if(method == 'asymmetric'){  # Kimber's boxplot
        vql <- qq[2] - qq[1]
        if(abs(vql) < 1e-06) warning("(Q2-Q1) = 0")
        vqu <- qq[3] - qq[2]
        if(abs(vqu) < 1e-06) warning("(Q3-Q2) = 0")
        low.b <- qq[1] - 2 * k * vql
        up.b  <- qq[3] + 2 * k * vqu
    }

    if(method == 'asymmetric2'){  # Walker et al. boxplot
        iqr <- qq[3] - qq[1]
        if(abs(iqr) < 1e-06){
            warning("IQR is 0: Bowley measure undefined, fences set at the quartiles")
            bow <- NA_real_
            vql <- vqu <- 0
        } else {
            bow <- (qq[3] + qq[1] - 2 * qq[2]) / iqr
            if(abs(bow) > 1 - 1e-06){
                warning("|Bowley skewness| = 1 (Q2 equals Q1 or Q3): ",
                        "symmetric IQR-based fences used instead")
                vql <- vqu <- iqr
            } else {
                vql <- iqr * (1 - bow) / (1 + bow)
                vqu <- iqr * (1 + bow) / (1 - bow)
            }
        }
        if(verbose) message('The Bowley skewness measure is: ', round(bow, 4))
        low.b <- qq[1] - k * vql
        up.b  <- qq[3] + k * vqu
    }

    if(method == 'adjbox'){  # Hubert & Vandervieren adjusted boxplot (unweighted)
        # medcouple: retry with doScale=TRUE if the default call fails
        warning("With method='adjbox' the argument k is set equal to 1.5")
        medc <- try(robustbase::mc(yy, doScale = FALSE), silent = TRUE)
        ck <- inherits(medc, "try-error")
        if(ck){
            medc <- try(robustbase::mc(yy, doScale = TRUE), silent = TRUE)
            if(inherits(medc, "try-error"))
                stop("robustbase::mc() failed: ",
                     conditionMessage(attr(medc, "condition")))
        }
        if(verbose) message('The MedCouple skewness measure is: ', round(medc, 4))

        aa <- robustbase::adjboxStats(x = yy, coef = 1.5, doScale = ck)
        low.b <- aa$fence[1]
        up.b  <- aa$fence[2]
    }

    # --- outliers ----------------------------------------------------------
    low.b <- unname(low.b)
    up.b  <- unname(up.b)
    fences <- c(lower = low.b, upper = up.b)

    # strict inequalities: values exactly on a fence are not flagged
    lower <- yy < low.b
    upper <- yy > up.b
    outl  <- lower | upper

    if(verbose){
        if(!any(outl)) message('No outliers found')
        else{
            message('No. of outliers in left tail: ', sum(lower))
            message('No. of outliers in right tail: ', sum(upper), '\n')
        }
    }

    # consistent output structure whether or not outliers exist
    list(quartiles = qq, fences = fences,
         excluded = to.check, outliers = lab[outl],
         lowOutl = lab[lower], upOutl = lab[upper])
}
