plot4ratios <- function(out){
    # plotting ratios=yt2/yt1 versus sizeU
    
    # check output
    if(is.null(out[["data"]])) stop("\n The input is not the output of HBmethod() or ratioSize(), or the functions were launched with the argument return.dataframe=FALSE")
    # get output data output of HBmethod() or ratioSize()
    dd <- out[["data"]]
    median.r <- out[["median.r"]]
    
    #check the type of out
    fn <- strsplit(x = as.character(out$call), split = "\\(")
    if(fn[[1]] == "HBmethod"){
        low.b <- out[["bounds.E"]][1]
        up.b <- out[["bounds.E"]][2]    
        rr <- dd$ratio
        # calculate bounds for sizeU
        Ls <- median.r * 1/(1 - low.b/dd$sizeU)
        Us <- median.r * (1 + up.b/dd$sizeU)
        lab <- "Ratios"
    }
    if(fn[[1]] ==  "ratioSize"){
        low.b <- out[["bounds"]][1]
        up.b <- out[["bounds"]][2]    
        # calculate bounds for sizeU
        rr <- dd$c.ratio
        Ls <- rep(low.b, nrow(dd))
        Us <- rep(up.b, nrow(dd))
        lab <- "centered Ratios"
    }
    # reorder units
    mat <- data.frame(id=dd$id, x=dd$sizeU, y=rr, yLow=Ls, yUp=Us, outl=dd$outliers)
    mat <- mat[order(mat$x), ]
    
    # draw scatterplot and bounds
    plot(x=mat$x, y=mat$y, type="n", 
         ylab = lab, xlab = "sizeU",main=paste(lab, "SizeU", sep=" Vs. "),
         ylim=c(min(c(mat$x, mat$yLow))-5, max(c(mat$x, mat$yUp))))
    
    
    tst <- mat$outl==1
    
    points(x=mat$x[tst],  y=mat$y[tst],  pch=18, col=2)
    points(x=mat$x[!tst], y=mat$y[!tst], pch=16, col=3)
    text(x=mat$x[tst], y=mat$y[tst], labels = mat$id[tst], 
         pos=1, cex=0.7, col=2)
    lines(x=mat$x, y=mat$yLow, col=4, lty=1, lwd=1.5)
    lines(x=mat$x, y=mat$yUp, col=4, lty=1, lwd=1.5)
    #legend("left", c("Outlier", "Not Outlier"), pch=c(18,16), col=c(2,3) )
    grid()
    colnames(mat) <- c("id", gsub(pattern = " ", replacement = ".", lab), "size","yLow", "yUp", "outlier")       
    mat
}