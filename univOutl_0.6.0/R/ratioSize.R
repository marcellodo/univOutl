ratioSize <- function(numerator, denominator, id = NULL,
                      size = NULL, U = 1, size.th = NULL,
                      return.dataframe = FALSE,
                      method = c("asymmetric2", "asymmetric", "adjbox", "resistant"),
                      ...) {
    
    method <- match.arg(method)
    
    # ---- input checks --------------------------------------------------------
    n <- length(numerator)
    if (n != length(denominator)) stop("Input vectors have different lengths")
    if (is.null(id)) id <- seq_len(n)
    if (length(id) != n) stop("'id' must have the same length as 'numerator'")
    if (!is.null(size) && length(size) != n)
        stop("'size' must have the same length as 'numerator'")
    
    # ---- discard NAs and zeros (NA-safe) --------------------------------------
    tst.na <- is.na(numerator) | is.na(denominator)
    tst.0  <- (numerator == 0) | (denominator == 0)
    tst.0[is.na(tst.0)] <- FALSE
    tst  <- tst.0 | tst.na
    keep <- !tst
    
    if (!any(keep)) stop("No valid units left after removing NAs and zeros")
    
    yy0 <- denominator[keep]
    yy1 <- numerator[keep]
    lab <- id[keep]
    
    if (any(yy0 < 0 | yy1 < 0))
        stop("Negative values found: the score is defined for positive ratios only")
    if (length(yy1) < 3)
        stop("At least 3 valid units are needed to identify outliers")
    
    # excluded units
    if (return.dataframe) {
        excluded <- data.frame(id = id[tst],
                               numerator = numerator[tst],
                               denominator = denominator[tst])
    } else {
        excluded <- id[tst]            # integer(0)-like when nothing is excluded
    }
    
    # ---- derive scores --------------------------------------------------------
    rr  <- yy1 / yy0
    mdn <- median(rr)
    sc  <- ifelse(rr < mdn, 1 - mdn / rr, rr / mdn - 1)
    
    # ---- boxplot-based fences --------------------------------------------------
    sk <- skew.misc(x = sc)
    message("Bowley skewness measure of centered ratios: ", round(sk[3], 4))
    
    aa <- boxB(x = sc, method = method, ...)
    
    bounds <- as.numeric(aa$fences[1:2])
    names(bounds) <- c("lower", "upper")
    if (anyNA(bounds) || any(is.infinite(bounds)))
        stop("Fences could not be computed (too many identical scores?)")
    
    outl <- (sc < bounds[["lower"]]) | (sc > bounds[["upper"]])
    
    # ---- importance measure (always computed, aligned with filtered data) ----
    imp <- if (is.null(size)) pmax(yy0, yy1)^U else size[keep]^U
    
    sub <- data.frame(id = lab,
                      numerator = yy1,
                      denominator = yy0,
                      ratio = rr,
                      c.ratio = sc,
                      outliers = as.integer(outl),
                      sizeU = imp)
    sub <- sub[order(sub$sizeU, decreasing = TRUE), ]
    
    # outliers reported to the user (optionally above the size threshold)
    sel <- sub$outliers == 1
    if (!is.null(size.th)) sel <- sel & (sub$sizeU > size.th^U)
    
    if (!any(outl)) message("No outliers found")
    
    # ---- output ----------------------------------------------------------------
    fine <- list(median.r = mdn, bounds = bounds, excluded = excluded)
    
    if (return.dataframe) {
        fine$data <- sub
    } else {
        fine$outliers <- sub$id[sel]
    }
    
    c(fine, list(call = match.call()))
}