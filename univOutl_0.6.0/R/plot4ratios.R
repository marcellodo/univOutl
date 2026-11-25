
plot4ratios <- function(out) {
    
    # Check input
    if (is.null(out[["data"]])) {
        stop("The input must be the output of HBmethod() or ratioSize(), ",
             "with return.dataframe = TRUE.")
    }
    
    dd <- out[["data"]]
    median.r <- out[["median.r"]]
    
    # Identify the originating function
    fn <- as.character(out$call)[1]
    fn <- sub("\\(.*$", "", fn)
    
    # Check required columns and calculate bounds
    if (fn == "HBmethod") {
        
        required <- c("id", "sizeU", "ratio", "outliers")
        missing <- setdiff(required, names(dd))
        
        if (length(missing) > 0)
            stop("Missing columns: ", paste(missing, collapse = ", "))
        
        low.b <- out[["bounds.E"]][1]
        up.b  <- out[["bounds.E"]][2]
        
        if (any(dd$sizeU <= low.b, na.rm = TRUE))
            warning("Some sizeU values are <= lower bound E.")
        
        rr <- dd$ratio
        Ls <- median.r / (1 - low.b/dd$sizeU)
        Us <- median.r * (1 + up.b/dd$sizeU)
        
        lab <- "Ratios"
        
    } else if (fn == "ratioSize") {
        
        required <- c("id", "sizeU", "c.ratio", "outliers")
        missing <- setdiff(required, names(dd))
        
        if (length(missing) > 0)
            stop("Missing columns: ", paste(missing, collapse = ", "))
        
        low.b <- out[["bounds"]][1]
        up.b  <- out[["bounds"]][2]
        
        rr <- dd$c.ratio
        Ls <- rep(low.b, nrow(dd))
        Us <- rep(up.b, nrow(dd))
        
        lab <- "centered Ratios"
        
    } else {
        stop("Unknown input: expected output from HBmethod() or ratioSize().")
    }
    
    # Check dimensions and numeric values
    if (length(rr) != nrow(dd) ||
        length(Ls) != nrow(dd) ||
        length(Us) != nrow(dd)) {
        stop("Inconsistent lengths in plotting data.")
    }
    
    if (!is.numeric(dd$sizeU))
        stop("sizeU must be numeric.")
    
    # Prepare plotting data
    mat <- data.frame(
        id = dd$id,
        x = dd$sizeU,
        y = rr,
        yLow = Ls,
        yUp = Us,
        outl = dd$outliers
    )
    
    mat <- mat[order(mat$x), ]
    
    # Calculate y-axis limits
    yr <- range(c(mat$y, mat$yLow, mat$yUp), finite = TRUE)
    
    if (!all(is.finite(yr)))
        stop("No finite values available for plotting.")
    
    if (diff(yr) == 0) {
        yr <- yr + c(-1, 1) * max(abs(yr), 1) * 0.05
    } else {
        yr <- yr + c(-1, 1) * diff(yr) * 0.05
    }
    
    # Draw scatterplot
    plot(
        x = mat$x,
        y = mat$y,
        type = "n",
        ylab = lab,
        xlab = "sizeU",
        main = paste(lab, "SizeU", sep = " Vs. "),
        ylim = yr
    )
    
    # Identify outliers
    tst <- !is.na(mat$outl) & mat$outl == 1
    not.tst <- !is.na(mat$outl) & mat$outl != 1
    
    # Plot observations
    points(mat$x[tst], mat$y[tst], pch = 18, col = 2)
    points(mat$x[not.tst], mat$y[not.tst], pch = 16, col = 3)
    
    # Label outliers
    text(
        mat$x[tst],
        mat$y[tst],
        labels = mat$id[tst],
        pos = 1,
        cex = 0.7,
        col = 2
    )
    
    # Plot bounds
    lines(mat$x, mat$yLow, col = 4, lty = 1, lwd = 1.5)
    lines(mat$x, mat$yUp, col = 4, lty = 1, lwd = 1.5)
    
    grid()
    
    # Return plotting data
    colnames(mat) <- c("id", "sizeU", lab, "yLow", "yUp", "outlier")
    
    return(mat)
}
