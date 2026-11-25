
HBmethod <- function(yt1, yt2, U = 0.5, A = 0.05, C = 4,
                     pct = 0.25, id = NULL, std.score = FALSE,
                     return.dataframe = FALSE) {
    
    # --------------------------------------------------
    # 1. Check input vectors
    # --------------------------------------------------
    
    if (!is.numeric(yt1) || !is.numeric(yt2))
        stop("yt1 and yt2 must be numeric vectors.")
    
    if (length(yt1) != length(yt2))
        stop("Input vectors have different lengths.")
    
    if (length(yt1) < 2)
        stop("At least two observations are required.")
    
    if (is.null(id)) {
        id <- seq_along(yt1)
    } else if (length(id) != length(yt1)) {
        stop("id must have the same length as yt1 and yt2.")
    }
    
    # --------------------------------------------------
    # 2. Check parameters
    # --------------------------------------------------
    
    if (length(U) != 1 || !is.numeric(U) ||
        !is.finite(U) || U <= 0 || U > 1)
        stop("U must be a numeric scalar strictly >0 and <=1.")
    
    if (length(A) != 1 || !is.numeric(A) ||
        !is.finite(A) || A <= 0)
        stop("A must be a positive finite numeric scalar.")
    
    if (length(pct) != 1 || !is.numeric(pct) ||
        !is.finite(pct) || pct <= 0 || pct >= 0.5)
        stop("pct must be a numeric scalar between 0 and 0.5.")
    
    if (length(C) == 1)
        C <- rep(C, 2)
    
    if (!is.numeric(C) || length(C) != 2 ||
        any(!is.finite(C)) || any(C < 0))
        stop("C must be a non-negative numeric scalar or vector of length 2.")
    
    if (!is.logical(std.score) || length(std.score) != 1 ||
        is.na(std.score))
        stop("std.score must be TRUE or FALSE.")
    
    if (!is.logical(return.dataframe) ||
        length(return.dataframe) != 1 ||
        is.na(return.dataframe))
        stop("return.dataframe must be TRUE or FALSE.")
    
    # if (!is.logical(adjboxE) || length(adjboxE) != 1 ||
    #     is.na(adjboxE))
    #     stop("adjboxE must be TRUE or FALSE.")
    
    # --------------------------------------------------
    # 3. Remove invalid observations
    # --------------------------------------------------
    
    tst <- !is.finite(yt1) | !is.finite(yt2) |
        yt1 <= 0 | yt2 <= 0
    
    discard <- data.frame(
        id = id[tst],
        yt1 = yt1[tst],
        yt2 = yt2[tst]
    )
    
    yy1 <- yt1[!tst]
    yy2 <- yt2[!tst]
    lab <- id[!tst]
    
    if (length(yy1) < 2)
        stop("Insufficient valid observations after exclusions.")
    
    # --------------------------------------------------
    # 4. Calculate Hidiroglou-Berthelot scores
    # --------------------------------------------------
    
    rr <- yy2 / yy1
    
    if (any(!is.finite(rr)) || any(rr <= 0))
        stop("Ratios must be positive and finite.")
    
    mdn.rr <- median(rr)
    
    if (!is.finite(mdn.rr) || mdn.rr <= 0)
        stop("Invalid median ratio.")
    
    sc <- ifelse(
        rr < mdn.rr,
        1 - mdn.rr / rr,
        rr / mdn.rr - 1
    )
    
    sizeU <- pmax(yy1, yy2)^U
    
    if (any(!is.finite(sizeU)))
        stop("Non-finite sizeU values. Check U and input data.")
    
    E <- sc * sizeU
    
    if (any(!is.finite(E)))
        stop("Non-finite E scores obtained.")
    
    # --------------------------------------------------
    # 5. Calculate quartiles of E scores
    # --------------------------------------------------
    
    q.E <- quantile(
        x = E,
        probs = c(pct, 0.50, 1 - pct),
        names = TRUE
    )
    
    if (all(abs(q.E) < 1e-6))
        stop("Quartiles of E scores are all equal to 0.")
    
    # --------------------------------------------------
    # 6. MedCouple & Bowley's skewness measure
    # --------------------------------------------------
    
    
    skew.E <- if (length(unique(E)) > 1 && length(E) >= 3) {
         skew.misc(x = E, weights = NULL)
    } else {
        NA_real_
    }
    mc.E <- skew.E[2]
    bow.E <- skew.E[3]
    
    message(
        "MedCouple skewness measure of E scores: ",
        round(mc.E, 4)
    )
    message(
        "Bowley's skewness measure of E scores: ",
        round(bow.E, 4)
    )
    
    # --------------------------------------------------
    # 7. Calculate asymmetric fences
    # --------------------------------------------------
    
    dq1 <- max(
        q.E[2] - q.E[1],
        abs(A * q.E[2])
    )
    
    dq3 <- max(
        q.E[3] - q.E[2],
        abs(A * q.E[2])
    )
    
    if (!is.finite(dq1) || !is.finite(dq3) ||
        dq1 <= 0 || dq3 <= 0)
        stop("Invalid scale parameters dq1 or dq3.")
    
    low.b <- q.E[2] - C[1] * dq1
    up.b  <- q.E[2] + C[2] * dq3
    
    names(low.b) <- "low"
    names(up.b) <- "up"
    
    bounds.E <- c(low.b, up.b)
    
    # --------------------------------------------------
    # 8. Standardized E scores (optional)
    # --------------------------------------------------
    
    if (std.score) {
        
        ncost <- qnorm(1 - pct)
        
        if (!is.finite(ncost) || ncost <= 0)
            stop("Invalid normalization constant.")
        
        std.Escore <- ifelse(
            E < q.E[2],
            (E - q.E[2]) / dq1 / ncost,
            (E - q.E[2]) / dq3 / ncost
        )
        
    }
    
    # --------------------------------------------------
    # 9. Identify outliers
    # --------------------------------------------------
    
    outl <- (E < low.b) | (E > up.b)
    
    if (sum(outl) == 0) {
        message("No outliers found.")
    } else {
        message(
            "Outliers found in the left tail: ",
            sum(E < low.b)
        )
        
        message(
            "Outliers found in the right tail: ",
            sum(E > up.b)
        )
    }
    
    # # --------------------------------------------------
    # # 10. Adjusted boxplot on E scores (optional)
    # # --------------------------------------------------
    # 
    # if (adjboxE) {
    #     
    #     out.bb <- boxB(
    #         x = E,
    #         k = 1.5,
    #         method = "adjbox"
    #     )
    #     
    #     # Identify outliers directly from the fences
    #     outlBB <- (E < out.bb$fences[1]) |
    #         (E > out.bb$fences[2])
    #     
    #     message(
    #         "Outliers found in the left tail with adj. boxplot on E: ",
    #         sum(E < out.bb$fences[1])
    #     )
    #     
    #     message(
    #         "Outliers found in the right tail with adj. boxplot on E: ",
    #         sum(E > out.bb$fences[2])
    #     )
    #     
    # }
    # 
    # --------------------------------------------------
    # 10. Prepare common output
    # --------------------------------------------------
    
    fine0 <- list(
        median.r = mdn.rr,
        quartiles.E = q.E,
        bounds.E = bounds.E,
        skew.E = c(mc.E, bow.E)
    )
    
    # --------------------------------------------------
    # 12. Prepare output with dataframe
    # --------------------------------------------------
    
    if (return.dataframe) {
        
        df.outl <- data.frame(
            id = lab,
            yt1 = yy1,
            yt2 = yy2,
            ratio = rr,
            sizeU = sizeU,
            Escore = E
        )
        
        if (std.score)
            df.outl$std.Escore <- std.Escore
        
        df.outl$outliers <- as.integer(outl)
        
        # if (adjboxE)
        #     df.outl$outliersBB <- as.integer(outlBB)
        
        fine1 <- list(
            excluded = discard,
            data = df.outl
        )
        
    } else {
        
        fine1 <- list(
            excluded = discard,
            outliers = lab[outl]
        )
        
        # if (adjboxE)
        #     fine1$outliersBB <- lab[outlBB]
        
    }
    
    # --------------------------------------------------
    # 13. Return results
    # --------------------------------------------------
    
    result <- c(fine0, fine1, list(call = match.call()))
    
    return(result)
}
