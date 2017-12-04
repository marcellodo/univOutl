HBmethod <- function(yt1, yt2, U=0.5, A=0.05, C=4, 
                    id=NULL, std.score=FALSE, return.dataframe=FALSE)
{

    # check input vectors
    if(length(yt1) != length(yt2)) stop('Input vectors have different lengths')
    if(is.null(id)) id <- 1:length(yt1)
                        
    # discard 0s and NAs
    tst.na <- is.na(yt1) | is.na(yt2)
    tst.0 <- (yt1==0) | (yt2==0)
    tst <- tst.0 | tst.na
    
    if(sum(tst)==0){
        discard <- integer(0)
        yy1 <- yt1
        yy2 <- yt2
        lab <- id
    }
    else{
        discard <- id[tst]
        yy1 <- yt1[!tst]
        yy2 <- yt2[!tst]
        lab <- id[!tst]
    }
                        
    # derive Hidiroglou-Berthelot score

    rr <- yy2/yy1
    mdn.rr <- median(rr)

    sc <- ifelse((rr < mdn.rr), (1 - mdn.rr/rr), (rr/mdn.rr - 1))
    E <- sc * (pmax(yy1, yy2)) ^ U
    
    # compute bounds for scores
    q.E <- quantile(x = E, probs = c(0.25, 0.50, 0.75))
    if(all(abs(q.E)<1e-06)) stop("Quartiles of E scores are all equal to 0")
    
    message('MedCouple skewness measure of E scores: ', round(robustbase::mc(E), 4))
    
    dq1 <- max( (q.E[2] - q.E[1]), abs(A * q.E[2]) )
    dq3 <- max( (q.E[3] - q.E[2]), abs(A * q.E[2]) )

    if(std.score){
        std.Escore <- (E - q.E[2])/dq1*0.6745
        std.Escore[E >= q.E[2]] <- c((E - q.E[2])/dq3*0.6745)[E >= q.E[2]]
    }

    # identifies outliers

    if(length(C)==1) C <- rep(C, 2)
    low.b <- q.E[2] - C[1] * dq1
    up.b <- q.E[2] + C[2] * dq3        

    names(low.b) <- 'low'
    names(up.b) <- 'up'
    
    outl <- (E < low.b) | (E > up.b)
    if(sum(outl)==0) message('No outliers found')
    else{
        message('Outliers found in the left tail: ', sum(E < low.b))
        message('Outliers found in the right tail: ', sum(E > up.b))
    }
    # output
    fine0 <- list(median.r = mdn.rr, quartiles.E=q.E, bounds.E=c(low.b, up.b))
    if(return.dataframe){
        discard <- data.frame(id=id[tst],
                              yt1=yt1[tst],
                              yt2=yt2[tst])

        df.outl <- data.frame(id=id[!tst],
                              yt1=yt1[!tst],
                              yt2=yt2[!tst],
                              ratio=rr, Escore=E)
        if(std.score) df.outl$std.Escore <- std.Escore
        df.outl$outliers <- as.integer(outl)
        
        fine1 <- list(excluded=discard, data = df.outl)
    }
    else{
        if(sum(outl)==0) fine1 <- list(excluded=discard, outliers = integer(0))   
        
        else fine1 <- list(excluded=discard, outliers = lab[outl])
    }
    c(fine0, fine1)
}