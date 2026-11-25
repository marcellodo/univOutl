
skew.misc <- function(x, weights = NULL) {
    
    # --------------------------------------------------
    # 1. Check input
    # --------------------------------------------------
    
    if (!is.numeric(x))
        stop("x must be a numeric vector.")
    
    if (length(x) < 2)
        stop("At least two observations are required.")
    
    # Remove missing and infinite observations
    valid <- is.finite(x)
    y <- x[valid]
    
    if (length(y) < 2)
        stop("Insufficient finite observations.")
    
    # --------------------------------------------------
    # 2. Check weights
    # --------------------------------------------------
    
    if (!is.null(weights)) {
        
        if (!is.numeric(weights) ||
            length(weights) != length(x))
            stop("weights must be numeric and have the same length as x.")
        
        if (any(!is.finite(weights)) || any(weights < 0))
            stop("weights must be finite and non-negative.")
        
        ww <- weights[valid]
        
        if (sum(ww) <= 0)
            stop("Sum of weights must be positive.")
        
        if (sum(ww > 0) < 2)
            stop("At least two observations must have positive weights.")
        
    }
    
    # --------------------------------------------------
    # 3. Calculate percentile-based skewness
    # --------------------------------------------------
    
    probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
    
    if (is.null(weights)) {
        
        qq <- quantile(
            x = y,
            probs = probs,
            names = FALSE
        )
        
    } else {
        
        qq <- Hmisc::wtd.quantile(
            x = y,
            weights = ww,
            probs = probs,
            normwt = FALSE,
            na.rm = TRUE
        )
        
        qq <- as.numeric(qq)
        
    }
    
    # Extract quantiles
    q10 <- qq[1]
    q25 <- qq[2]
    q50 <- qq[3]
    q75 <- qq[4]
    q90 <- qq[5]
    
    # Bowley quartile skewness
    den.Q <- q75 - q25
    
    skQ <- if (den.Q > 0) {
        (q75 - 2 * q50 + q25) / den.Q
    } else {
        NA_real_
    }
    
    # Decile-based skewness
    den.P <- q90 - q10
    
    skD <- if (den.P > 0) {
        (q90 - 2 * q50 + q10) / den.P
    } else {
        NA_real_
    }
    
    # Transformations
    gQ <- if (is.finite(skQ) && abs(1 - skQ) > 1e-10) {
        (1 + skQ) / (1 - skQ)
    } else {
        NA_real_
    }
    
    gD <- if (is.finite(skD) && abs(1 - skD) > 1e-10) {
        (1 + skD) / (1 - skD)
    } else {
        NA_real_
    }
    
    out <- c(
        Bowley.Q = skQ,
        g.Q = gQ,
        Bowley.P = skD,
        g.P = gD
    )
    
    # --------------------------------------------------
    # 4. MedCouple
    # --------------------------------------------------
    
    if (length(y) >= 3 && length(unique(y)) > 1) {
        
        medc <- tryCatch(
            robustbase::mc(y, doScale = FALSE),
            error = function(e) {
                tryCatch(
                    robustbase::mc(y, doScale = TRUE),
                    error = function(e2) NA_real_
                )
            }
        )
        
    } else {
        medc <- NA_real_
    }
    
    # --------------------------------------------------
    # 5. Pearson moment coefficient
    # --------------------------------------------------
    
    z <- y - mean(y)
    
    # Rescale to improve numerical stability
    zmax <- max(abs(z))
    
    if (zmax > 0) {
        
        z <- z / zmax
        
        pears <- mean(z^3) /
            (mean(z^2)^(3/2))
        
    } else {
        pears <- NA_real_
    }
    
    # --------------------------------------------------
    # 6. Return results
    # --------------------------------------------------
    
    c(
        Pearson = pears,
        MedCouple = as.numeric(medc),
        out
    )
}
