ratioSize <- function(numerator, denominator, id=NULL,
                      size=NULL, U=1, size.th=NULL, return.dataframe=FALSE){
    
    # check lengt of vectors and identifiers
    if(length(numerator) != length(denominator)) stop("Input vectors have different lengths")
    if(is.null(id)) id <- 1:length(numerator)
    
    # discard 0s and NAs
    tst.na <- is.na(numerator) | is.na(denominator)
    tst.0 <- (numerator==0) | (denominator==0)
    tst <- tst.0 | tst.na
    #message('There are ', sum(tst), ' units with NAs or 0s')
    
    if(sum(tst)==0){
        df.no.outl <- integer(0)
        yy0 <- denominator
        yy1 <- numerator
        lab <- id
    }
    else{
        df.no.outl <- id[tst]
        yy0 <- denominator[!tst]
        yy1 <- numerator[!tst]
        lab <- id[!tst]
    }
    
    # derive scores
    rr <- yy1/yy0
    mdn.ratios <- median(rr)
   
    sc <- ifelse((rr < mdn.ratios), (1 - mdn.ratios/rr), (rr/mdn.ratios - 1))
    
    # computes ranges using adjusted boxplot
    message('MedCouple skewness measure of centerad ratios: ', round(robustbase::mc(sc), 4))
    aa <- robustbase::adjboxStats(x = sc)
    
    # identifies outliers
    outl <- (sc < aa$fence[1]) | (sc > aa$fence[2])
    
    sub <- data.frame(id=lab,
                      numerator=yy1,
                      denominator=yy0,
                      ratio=rr, c.ratio=sc, 
                      outliers=as.integer(outl))
    
    if(sum(outl)==0) message('No outliers found')
    # else{
    #     message('No. of outliers in left tail: ', sum(sc < aa$fence[1]))
    #     message('No. of outliers in right tail: ', sum(sc > aa$fence[2]), '\n')
    # }
    
    # output 
    
    fine0 <- list(median.r=mdn.ratios, bounds = c(aa$fence))
    
    if(return.dataframe){
        df.no.outl <- data.frame(id=id[tst],
                              numerator=numerator[tst],
                              denominator=denominator[tst])
    }
    # when no outliers are found
    if(sum(outl)==0){
        if(return.dataframe){
            fine <- c(fine0, list(excluded=df.no.outl, data=sub))
        }
        else{
            fine <- c(fine0, list(excluded=df.no.outl, 
                         outliers=integer(0)))    
        } 

    }
    # output when outliers are found
    # importance measure 'size' is considered 
    # when available also a threshold 'size.th'
    else{
        if(is.null(size) ) imp <- pmax(yy0, yy1) ^ U
        else imp <- size^U
        sub$sizeU <- imp
        sub <- sub[order(imp, decreasing=TRUE), ]
        # sub <- sub[sub$outlier, ]
        
        if(!is.null(size.th)) {
            sub <- sub[((sub$sizeU > (size.th^U)) & (sub$outliers==1)), ]
        }
        if(return.dataframe){
            fine <- c(fine0, list(excluded=df.no.outl, 
                         data=sub))
        }
        else {
            fine <- c(fine0, list(excluded=df.no.outl, 
                     outliers=sub$id[sub$outliers==1]))
        }
    }
    c(fine, call=match.call())
}