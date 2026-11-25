\name{boxB}
\alias{boxB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	BoxPlot based outlier detection
}
\description{
	Identifies univariate outliers by using methods based on BoxPlots
}
\usage{
boxB(x, k=1.5, method=c('asymmetric', 'asymmetric2', 'resistant', 
     'standard', 'adjbox'), weights=NULL, id=NULL, exclude=NA, 
     verbose=TRUE)
}

\arguments{
  \item{x}{
		Numeric data vector on which outlier detection will be performed. 
}
  \item{k}{
	Single nonnegative constant that determines the extension of the 'whiskers'. Commonly used values are 1.5 (default, inner fence), or 3 (outer fence, rarely used). The constant \code{k} is used by all the methods;but, with \code{method="adjbox"} only the value 1.5 is allowed (see Details).
}
  \item{method}{
		Character, identifies the method to be used: \code{method="resistant"} (or equivalently \code{method="standard"}) provides the `standard' boxplot fences; \code{method="asymmetric"} is a modification of standard method (Kimber, 1990) to deal with (moderately) skewed data; \code{method="asymmetric2"} is the modification proposed by Walker et al (2018) to deal with asymmetric distributions;\code{method="adjbox"} uses Hubert and Vandervieren (2008) adjusted boxplot for skewed distributions. 
}
  \item{weights}{
		Optional numeric vector with units' weights associated to the observations in \code{x}; it must have the same length of \code{x}. Only nonnegative weights (without missing values) are allowed. Weights are used in estimating the quartiles (see Details). Note that they are ignored (with a warning) when \code{method="adjbox"}.
}
  \item{id}{
Optional vector with identifiers of units in \code{x}; it must have the same length of \code{x}. If missing (\code{id=NULL}, default) the identifiers will be set equal to the positions in the vector (i.e. \code{id=1:length(x)}). 
}
\item{exclude}{
Values of \code{x} that will be excluded by the analysis. Missing values (\code{NA} and \code{NaN}) are always excluded, whatever the value of this argument (default is \code{exclude=NA}). 
}
  \item{verbose}{
Logical, if \code{TRUE} (default) messages are printed (skewness measure, number of outliers in each tail, etc.). Set it to \code{FALSE} to suppress them, e.g. when the function is called in a loop.
}

}
\details{
When \code{method="resistant"} (or \code{"standard"}) the outlying observations are those outside the interval:

\deqn{[Q_1 - k \times IQR;\quad  Q_3 + k \times IQR] }{[Q1 - k*IQR; Q3 + k*IQR]}

where \eqn{Q_1}{Q1} and \eqn{Q_3}{Q3} are respectively the 1st and the 3rd quartile of \code{x}, while \eqn{IQR=(Q_3 - Q_1)}{IQR=(Q3-Q1)} is the Inter-Quartile Range. The value \eqn{k=1.5}{k=1.5} (said 'inner fences') is commonly used when drawing a boxplot. Values \eqn{k=2}{k=2} and  \eqn{k=3}{k=3} provide middle and outer fences, respectively.

When \code{method="asymmetric"} the outlying observations are those outside the interval: 

\deqn{[Q_1 - 2k \times (Q_2-Q_1);\quad  Q_3 + 2k \times (Q_3-Q_2)] }{[Q1 - 2k*(Q2-Q1); Q3 + 2k*(Q3-Q2)]}

being \eqn{Q_2}{Q2} the median; such a modification allows to account for slight skewness of the distribution. A warning is issued when \eqn{Q_2-Q_1}{Q2-Q1} or \eqn{Q_3-Q_2}{Q3-Q2} is (almost) zero, because in such a case the corresponding fence collapses on the quartile.

When \code{method="asymmetric2"} the outlying observations are identified using the Walker et al (2018) modified boxplot: 

\deqn{[Q_1 - k \times IQR \times \frac{(1-B)}{(1+B)};\quad  Q_3 + k \times IQR \times \frac{(1+B)}{(1-B)}] }{[Q1 - k*IQR*(1-B)/(1+B); Q3 + k*IQR*(1+B)/(1-B)]}

being B the Bowley's skewness coefficient (see \code{\link[univOutl]{skew.misc}}). In the degenerate cases the function issues a warning: when \eqn{IQR=0}{IQR=0} the Bowley's coefficient is undefined and the fences are set equal to the quartiles; when \eqn{|B|=1}{|B|=1} (i.e. \eqn{Q_2}{Q2} equals \eqn{Q_1}{Q1} or \eqn{Q_3}{Q3}) the symmetric fences \eqn{Q_1 - k \times IQR}{Q1 - k*IQR} and \eqn{Q_3 + k \times IQR}{Q3 + k*IQR} are used.

Finally, when \code{method="adjbox"} the outlying observations are identified using the method proposed by Hubert and Vandervieren (2008) and based on the Medcouple measure of skewness; in practice the bounds are: 

\deqn{[Q_1-k \times e^{aM} \times IQR;\quad  Q_3+k \times e^{bM}\times IQR ]}{[Q1 - k*exp(a*M)*IQR; Q3 + k*exp(b*M)*IQR]}

Where M is the medcouple; when \eqn{M > 0}{M>0} (positive skewness) then \eqn{a = -4}{a=-4} and \eqn{b = 3}{b=3}; on the contrary \eqn{a = -3}{a = -3} and \eqn{b = 4}{b = 4} for negative skewness (\eqn{M < 0}{M < 0}). This adjustment of the boxplot, according to Hubert and Vandervieren (2008), works with moderate skewness (\eqn{-0.6 \leq M \leq 0.6}{-0.6 <= M <= 0.6}). The bounds of the adjusted boxplot are derived by applying the function \code{\link[robustbase]{adjboxStats}} in the package \pkg{robustbase}, with \code{coef=1.5}. The medcouple is computed by \code{\link[robustbase]{mc}}; if the default call fails, the function tries again with \code{doScale=TRUE}.

When weights are available (passed via the argument \code{weights}) then they are used in the computation of the quartiles for the methods \code{"resistant"}, \code{"standard"}, \code{"asymmetric"} and \code{"asymmetric2"}. In particular, the quartiles are derived using the function  \code{\link[Hmisc]{wtd.quantile}} in the package \pkg{Hmisc}, with weights treated as frequencies (\code{normwt=FALSE}). Note that the definition of the (weighted) quantiles used by \code{\link[Hmisc]{wtd.quantile}} differs slightly from the one used by \code{\link[stats]{quantile}} (type 7) in the unweighted case, so the results may differ even when all the weights are equal.

The medcouple has no weighted version, therefore with \code{method="adjbox"} the weights are ignored (a warning is issued).


}

\value{
The output is a list containing the following components:

\item{quartiles}{
The quartiles of \code{x} after discarding the values in the \code{exclude} argument (and the missing values). When weights are provided (and \code{method} is not \code{"adjbox"}) they are used in quartiles estimation through the function  \code{\link[Hmisc]{wtd.quantile}} in the package \pkg{Hmisc}.
}

\item{fences}{
The bounds of the interval (named \code{lower} and \code{upper}), values outside the interval are detected as outliers.
}

\item{excluded}{
	The identifiers or positions (when \code{id=NULL}) of units in \code{x} excluded by the computations, according to the argument \code{exclude} (missing values included).
}

\item{outliers}{
The identifiers or positions (when \code{id=NULL}) of units in \code{x} detected as outliers.
}

\item{lowOutl}{
The identifiers or positions (when \code{id=NULL}) of units in \code{x} detected as outliers in the lower tail of the distribution.
}

\item{upOutl}{
The identifiers or positions (when \code{id=NULL}) of units in \code{x} detected as outliers in the upper tail of the distribution.
}

All the components are always returned; when no outliers are found \code{outliers}, \code{lowOutl} and \code{upOutl} are empty vectors.

}
\references{
McGill, R., Tukey, J. W. and Larsen, W. A. (1978) `Variations of box plots'. \emph{The American Statistician}, 32, pp. 12-16.

Kimber, A.C. (1990) `Exploratory data analysis for possibly censored data from skewed distributions'. \emph{Applied Statistics}, 39, pp. 21-30.

Hubert, M., and Vandervieren, E. (2008) `An Adjusted Boxplot for Skewed Distributions', \emph{Computational Statistics and Data Analysis}, 52, pp. 5186-5201.

Walker, M.L., Dovoedo, Y.H., Chakraborti, S., and Hilton, C.W. (2018) `An Improved Boxplot for Univariate Data'. \emph{The American Statistician}, 72, pp. 348-353

}
\author{
Marcello D'Orazio \email{mdo.statmatch@gmail.com}
}


\seealso{
\code{\link[robustbase]{adjboxStats}}, \code{\link[robustbase]{mc}}, \code{\link[Hmisc]{wtd.quantile}}, \code{\link[univOutl]{skew.misc}} 
}
\examples{

set.seed(321)
x <- rnorm(30, 50, 10)
x[10] <- 1
x[20] <- 100

out <- boxB(x = x, k = 1.5, method = 'asymmetric')
out$fences
out$outliers
x[out$outliers]

out <- boxB(x = x, k = 1.5, method = 'asymmetric2')
out$fences
out$outliers
x[out$outliers]

out <- boxB(x = x, k = 1.5, method = 'adjbox')
out$fences
out$outliers
x[out$outliers]

# missing values are excluded
x[24] <- NA
x.ids <- paste0('obs',1:30)
out <- boxB(x = x, k = 1.5, method = 'adjbox', id = x.ids)
out$excluded
out$fences
out$outliers

# weights are used in the estimation of the quartiles
set.seed(111)
w <- round(runif(n = 30, min=1, max=10))
out <- boxB(x = x, k = 1.5, method = 'asymmetric', id = x.ids, weights = w)
out$excluded
out$quartiles
out$fences
out$outliers
out$lowOutl
out$upOutl

# with method='adjbox' weights are ignored (a warning is issued)
out <- boxB(x = x, k = 1.5, method = 'adjbox', id = x.ids, weights = w, 
            verbose = FALSE)
out$fences

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ survey }
