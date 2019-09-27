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
        lab <- "cent. ratios"
    }
    
    
    # draw scatterplot and bounds
    plot(x=dd$sizeU, y=rr, type="n", 
         ylab = lab, xlab = "sizeU",main=paste(lab, "SizeU", sep=" Vs. "),
         ylim=c(min(c(dd$ratio, Ls)),max(c(dd$ratio,Us))))
    
    lines(x=dd$sizeU, y=Ls, col=4, lty=1, lwd=1.5)
    lines(x=dd$sizeU, y=Us, col=4, lty=1, lwd=1.5)
    
    tst <- dd$outliers==1
    
    points(x=dd$sizeU[tst], y=rr[tst], pch=18, col=2, cex=1.5)
    points(x=dd$sizeU[!tst], y=rr[!tst], pch=16, col=3)
    text(x=dd$sizeU[tst], y=rr[tst], labels = dd$id[tst], 
         pos=1, cex=1, col=2)
    #legend("left", c("Outlier", "Not Outlier"), pch=c(18,16), col=c(2,3) )
    grid()
    list(x=dd$sizeU, y=dd$ratio, y.low=Ls, y.up=Us)       
    
}