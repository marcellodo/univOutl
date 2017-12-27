LocScaleB <- function(x, k=3, method='MAD',  weights=NULL, id=NULL, 
                      exclude=NA, logt=FALSE, return.dataframe=FALSE){
    
    # units' identifiers
    if(is.null(id)) id <- 1:length(x)
    
    # removes values in exclude vector
    tst <- x %in% exclude
    if(sum(tst)==0){
        to.check <- integer(0)
        yy <- x
        lab <- id
        ww <- weights
    }
    else{
        to.check <- id[tst]
        yy <- x[!tst]
        lab <- id[!tst]    
        ww <- weights[!tst]
    }

    # log transform variable
    
    if(logt){ 
        yy <- log(yy+1)
        warning('Please note that log(x+1) transformation is applied')
    }    
    
    # computes quantiles
    if(is.null(weights)) qq <- quantile(x=yy, probs=c(0.25, 0.50, 0.75))
    else qq <- Hmisc::wtd.quantile(x=yy, weights=ww, probs=c(0.25, 0.50, 0.75))
    
    if(all(abs(qq)<1e-06)) stop("Quartiles are all equal to 0")
    
    # estimates scale measure
    if(tolower(method) == 'idr'){
        if(is.null(weights)) qq10 <- quantile(x=yy, probs=c(0.10, 0.90))
        else qq10 <- Hmisc::wtd.quantile(x=yy, weights=ww, probs=c(0.10, 0.90))
        ddl <- ddr <- (qq10[2] - qq10[1])/2.5631
    }
    if(tolower(method) == 'iqr'){
        ddl <- ddr <- (qq[3] - qq[1])/1.3490
    }
    if(tolower(method) == 'dq'){
        ddl <- (qq[2] - qq[1])/0.6745
        ddr <- (qq[3] - qq[2])/0.6745 
    }
    if(tolower(method) == 'adjout'){
        medc <- robustbase::mc(yy)
        message('The MedCouple skewness measure is: ', round(medc,4))
        if(is.null(weights)){
            aa <- robustbase::adjboxStats(x=yy)
            low <- aa$fence[1]
            up <- aa$fence[2]    
        }
        else{
            if(medc >= 0){
                low <- qq[1] - 1.5*exp(-4*medc)*(qq[3]-qq[1])
                up  <- qq[3] + 1.5*exp( 3*medc)*(qq[3]-qq[1])
            }
            else{
                low <- qq[1] - 1.5*exp(-3*medc)*(qq[3]-qq[1])
                up  <- qq[3] + 1.5*exp( 4*medc)*(qq[3]-qq[1])
            }
        }
        ddl <- qq[2] - low
        ddr <- up - qq[2]
    }
    if(tolower(method) == 'mad'){
        if(is.null(weights)){
            ddl <- ddr <- mad(yy)
        }
        else{
            zz <- abs(yy - qq[2])
            ddl <- ddr <- Hmisc::wtd.quantile(x=zz, weights=ww, probs=0.5) * 1.4826
        }
    }
    if(tolower(method) == 'gini'){
        if(is.null(weights)){
            ddl <- ddr <- Hmisc::GiniMd(yy) / 2 * sqrt(pi)
        }
        else{
            stop('Weights cannot be used when computing Gini Mean Difference')
        }
    }
    if(tolower(method) == 'scaletau2'){
        if(is.null(weights)){
            ddl <- ddr <- robustbase::scaleTau2(yy)
        }
        else{
            stop('Weights cannot be used when computing scaletau2')
        }
    }
    if(tolower(method) == 'qn'){
        if(is.null(weights)){
            ddl <- ddr <- robustbase::Qn(yy)
        }
        else{
            stop('Weights cannot be used when computing Qn')
        }
    }
    if(tolower(method) == 'sn'){
        if(is.null(weights)){
            ddl <- ddr <- robustbase::Sn(yy)
        }
        else{
            stop('Weights cannot be used when computing Sn')
        }
    }
    
    # check ddl and ddr
    if(ddl<=0){
        warning("Estimated scale measure is 0 \n
                it will be substituted with |0.05*Median|")
        ddl <- abs(0.05 * qq[2])
    }
    if(ddr<=0){
        warning("Estimated scale measure is 0 \n
                it will be substituted with |0.05*Median|")
        ddr <- abs(0.05 * qq[2])
    }
    # derive score
    ee <- yy - qq[2]
    if(ddl == ddr) zz <- ee/ddl
    else zz <- ifelse(yy < qq[2], ee/ddl, ee/ddr)
    
    # derives bounds
    low.b <- qq[2] - k*ddl 
    up.b <- qq[2] + k*ddr
    names(low.b) <- 'low'
    names(up.b) <- 'up'
    
    # identifies outliers
    outl <- (yy < low.b) | (yy > up.b)
    # outl <- (zz < -k) | (zz > k)
    if(sum(outl)==0) message('No outliers found')
    else{
        # message('No. of outliers in left tail: ', sum(zz < -k))
        # message('No. of outliers in right tail: ', sum(zz > k), '\n')
        message('No. of outliers in left tail: ', sum(yy < low.b))
        message('No. of outliers in right tail: ', sum(yy > up.b), '\n')
        
    }
        
    # output
    
    if(ddl==ddr){ 
        pp <- c(scale=unname(ddr))
    }
    else{
        pp <- c(sc.left=unname(ddl), sc.right=unname(ddr))  
    }
    pp <- c(median=c(unname(qq[2])), pp)
    
    fine0 <- list(pars=pp, bounds=c(lower=low.b, upper=up.b)) 
    
    if(return.dataframe){
        to.check <- data.frame(id=id[tst], x=x[tst])
        if(!is.null(weights)) to.check$weight <- weights[tst]
        
        df.out <- data.frame(id=lab, x=x[!tst])
        if(logt) df.out$log.x <- yy
        if(!is.null(weights)) df.out$weight <- ww
        df.out$score <- zz
        df.out$outliers <- as.integer(outl)
        fine1 <- list(excluded=to.check, data=df.out) 
    }
    else{
        if(sum(outl)==0){
            fine1 <- list(excluded=to.check, outliers=integer(0)) 
        }
        else{
            fine1 <- list(excluded=to.check, outliers=lab[outl]) 
        }    
    } 
    c(fine0, fine1)
}    