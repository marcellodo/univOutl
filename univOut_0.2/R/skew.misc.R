skew.misc <- function(x, weights=NULL){
    y <- x[!is.na(x)]
    if(!is.null(weights)) ww <- weights[!is.na(x)]
    
    # percentile-based skewness coefficient (Bowley)
    if(is.null(weights)) qq <- quantile(x = y, probs = c(0.10, 0.25, 0.50, 0.75, 0.90)) 
    else qq <- Hmisc::wtd.quantile(x=y, weights=ww, probs=c(0.10, 0.25, 0.50, 0.75, 0.90))
    
     
    skQ <- (qq["75%"] - 2*qq["50%"] + qq["25%"])/(qq["75%"] - qq["25%"])
    gQ<- (1+skQ)/(1-skQ)
    
    skD <- (qq["90%"] - 2*qq["50%"] + qq["10%"])/(qq["90%"] - qq["10%"])
    gD <- (1+skD)/(1-skD)
    
    out <- c(skQ, gQ, skD, gD)
    names(out) <- c("Bowley.Q", "g.Q", "Bowley.P", "g.P")
    
    # medcouple
    medc <- robustbase::mc(y)
    
    #pearson 
    m <- length(y)
    z <- y - mean(y)
    pears <- mean(z^3)/((mean(z^2))^(3/2))

    # output
    c(Pearson=pears, MedCouple=medc, out)
}