\name{ratioSize}
\alias{ratioSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifies outliers on ratios and filter them by a size measure
}
\description{
Identifies outliers on transformed ratios (centering with respect to their median) using the adjusted boxplot for skewed distributions. Outliers can be sorted/filtered according to a size measure.  
}
\usage{
ratioSize(numerator, denominator, id=NULL,
          size=NULL, U=1, size.th=NULL, return.dataframe=FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numerator}{
 Numeric vector with the values that go at numerator of the ratio
}
  \item{denominator}{
 Numeric vector with the values that go at denominator of the ratio
}
  \item{id}{
Optional numeric or character vector, with identifiers of units. If \code{id=NULL} units identifiers will be set equal to their positions in \code{x}.
}
  \item{size}{
Optional numeric vector providing a measure of the importance of a ratio. If \code{size = NULL} the size measure is the maximum value between the numerator and the denominator of each ratio (makes sense if both the variables are observed using the same unit of measure). Observations' importance is also controlled by the argument \code{U}.
}
  \item{U}{
Numeric, constant with \eqn{0 < U \leq 1}{0<U<=1} controlling importance of each unit, in practice the final size measure is derived as (size^U). Commonly used values are 0.4, 0.5 or 1.
}

  \item{size.th}{
Numeric, size threshold. Can be specified when a size measure is used. In such a case just outliers with a size greater than the threshold will be returned. Note that when argument \code{U} is not set equal to 1, then the final threshold will be \code{size.th^U}. 
}
 \item{return.dataframe}{
Logical, if \code{TRUE} the output will save all the relevant information for outlier detection in a dataframe with the following columns: `id' (units' identifiers), `numerator', `denominator', `ratio' (= numerator/denominator), `c.ratio' (centered ratios, see Details), `size' (size^U values) and finally `outliers', where value 1 indicates observations detected as an outlier and 0 otherwise.  
}


}
\details{

This function searches for outliers starting from ratios \code{r=numerator/denominator}. At first the ratios are centered around their median, as in Hidiroglou Berthelot (1986) procedure (see \code{\link[univOutl]{HBmethod}}), then the outlier identification is based on the adjusted boxplot for skewed distribution (Hubert and Vandervieren 2008) (see \code{\link[robustbase]{adjboxStats}}).
The subset of outliers is sorted in decreasing order according the size measure. If a size threshold is provided then just outliers with  (size^U) > (size.th^U)  will be returned. 

}
\value{

A list whose components depend on the \code{return.dataframe} argument. When \code{return.dataframe = FALSE} just the following components are returned:

\item{median.r}{the median of the ratios}
\item{bounds}{The bounds of the interval for centered ratios}

\item{excluded}{
	The position or the identifiers of the units with values excluded by the computations because of 0s or NAs.
}
\item{outliers}{
	The position or the identifiers of the units detected as outliers. Remember that when \code{size.th} is set, just outliers with  (size^U) > (size.th^U)  will be returned.  
}

When \code{return.dataframe=TRUE} the latter two components are substituted with two dataframes:

\item{excluded}{
	A dataframe with the subset of observations excluded 
}
	
\item{data}{
	A dataframe with the not excluded observations with the following columns: `id' (units' identifiers), `numerator', `denominator', `ratio' (= numerator/denominator), `c.ratio' (centered ratios, see Details), `size' (size^U values) and finally `outliers', where value 1 indicates observations detected as an outlier and 0 otherwise. The data frame will be sorted in decreasing manner according to size^U. Note that when a size threshold is provided then ONLY outliers with  (size^U) > (size.th^U)  will be returned. 
}


}
\references{
Hidiroglou, M.A. and Berthelot, J.-M. (1986) `Statistical editing and Imputation for Periodic Business Surveys'. \emph{Survey Methodology}, Vol 12, pp. 73-83.

Vanderviere, E. and Huber, M. (2008) `An Adjusted Boxplot for Skewed Distributions', \emph{Computational Statistics and Data Analysis}, 52, pp. 5186-5201.

}
\author{
Marcello D'Orazio \email{mdo.statmatch@gmail.com}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[univOutl]{HBmethod}}, \code{\link[univOutl]{boxB}},\code{\link[robustbase]{adjboxStats}}
}
\examples{

set.seed(444)
x1 <- rnorm(30, 50, 5)
set.seed(555)
rr <- runif(30, 0.9, 1.2)
rr[10] <- 2
x2 <- x1 * rr

out <- ratioSize(numerator = x2, denominator = x1)
out

out <- ratioSize(numerator = x2, denominator = x1,
                 return.dataframe = TRUE)
head(out$data)

out <- ratioSize(numerator = x2, denominator = x1,
                 size.th = 65, return.dataframe = TRUE)
head(out$data)


	
}
\keyword{ univar }
\keyword{ survey }
