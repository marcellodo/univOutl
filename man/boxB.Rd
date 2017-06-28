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
boxB(x, k=1.5, method='asymmetric', weights=NULL, id=NULL, 
     exclude=NA, logt=FALSE)
}

\arguments{
  \item{x}{
		Numeric vector that will be searched for outliers. 
}
  \item{k}{
	Nonnegative constant that determines the extension of the 'whiskers'. Commonly used values are 1.5 (default), 2, or 3.
Note that when \code{method="adjbox"} then \code{k} is set automatically equal to 1.5
}
  \item{method}{
		Character, identifies the method to be used: \code{method="resistant"} provides the `standard' boxplot fences; \code{method="asymmetric"} is a modification of standard method to deal with (moderately) skewed data; \code{method="adjbox"} uses Hubert and Vandervieren (2008) adjusted boxplot for skewed distributions.
}
  \item{weights}{
		Optional numeric vector with units' weights associated to observations of \code{x}. Only nonnegative weights are allowed. Weights are used in estimating the quartiles (see Details).
}
  \item{id}{
Optional vector with identifiers of units in \code{x}. If missing (\code{id=NULL}, default) the identifiers will be set equal to the positions in the vector (i.e. \code{id=1:length(x)}). 
}
\item{exclude}{
Values of \code{x} that will be excluded by the analysis. By default missing values are excluded (\code{exclude=NA}). 
}
  \item{logt}{
Logical, if \code{TRUE}, before searching outliers the \code{x} variable is log-transformed (\code{log(x+1)} is considered). In this case the summary outputs (bounds, etc.) will refer to the log-transformed \code{x}
}

}
\details{
When \code{method="resistant"} the outlying observations are those outside the interval:

\deqn{[Q_1 - k \times IQR;\quad  Q_3 + k \times IQR] }{[Q1 - k*IQR; Q3 + k*IQR]}

where \eqn{Q_1}{Q1} and \eqn{Q_3}{Q3} are respectively the 1st and the 3rd quartile of \code{x}, while \eqn{IQR=(Q_3 - Q_1)}{IQR=(Q3-Q1)} is the Inter-Quartile Range. The value \eqn{k=1.5}{k=1.5} (said 'inner fences') is commonly used when drawing a boxplot. Values \eqn{k=2}{k=2} and  \eqn{k=3}{k=3} provide middle and outer fences, respectively.

When \code{method="asymmetric"} the outlying observations are those outside the interval: 

\deqn{[Q_1 - 2k \times (Q_2-Q_1);\quad  Q_3 + 2k \times (Q_3-Q_2)] }{[Q1 - 2k*(Q2-Q1); Q3 + 2k*(Q3-Q2)]}

being \eqn{Q_2}{Q2} the median; such a modification allows to account for slight skewness of the distribution. 

Finally, when \code{method="adjbox"} the outlying observations are identified using the method proposed by Hubert and Vandervieren (2008) and based on the Medcouple measure of skewness; in practice the bounds are: 

\deqn{[Q_1-1.5 \times e^{aM} \times IQR;\quad  Q_3+1.5 \times e^{bM}\times IQR ]}{[Q1 - 1.5*exp(a*M)*IQR; Q3 + 1.5*exp(b*M)*IQR]}

Where M is the medcouple; when \eqn{M > 0}{M>0} (positive skewness) then \eqn{a = -4}{a=-4} and \eqn{b = 3}{b=3}; on the contrary \eqn{a = -3}{a = -3} and \eqn{b = 4}{b = 4} for negative skewness (\eqn{M < 0}{M < 0}). This adjustment of the boxplot, according to Hubert and Vandervieren (2008), works with moderate skewness (\eqn{-0.6 \leq M \leq 0.6}{-0.6 <= M <= 0.6}). The bounds of the adjusted boxplot are derived by applying the function \code{\link[robustbase]{adjboxStats}} in the package \pkg{robustbase}.

When weights are available (passed via the argument \code{weights}) then they are used in the computation of the quartiles. In particular, the quartiles are derived using the function  \code{\link[Hmisc]{wtd.quantile}} in the package \pkg{Hmisc}.  

Remember that when asking a log transformation (argument \code{logt=TRUE}) all the estimates (quartiles, etc.) will refer to \eqn{log(x+1)}{log(x+1)}.

}

\value{
The output is a list containing the following components:

\item{quartiles}{
The quartiles of \code{x} after discarding the values in the \code{exclude} argument. When weights are provided they are used in quartiles estimation trough the function  \code{\link[Hmisc]{wtd.quantile}} in the package \pkg{Hmisc}.
}

\item{fences}{
The bounds of the interval, values outside the interval are detected as outliers.
}

\item{excluded}{
	The identifiers or positions (when \code{id=NULL}) of units in \code{x} excluded by the computations, according to the argument \code{exclude}.
}

\item{outliers}{
The identifiers or positions (when \code{id=NULL}) of units in \code{x} detected as outliers.
}


}
\references{
McGill, R., Tukey, J. W. and Larsen, W. A. (1978) `Variations of box plots'. \emph{The American Statistician}, 32, pp. 12-16.

Vanderviere, E. and Huber, M. (2008) `An Adjusted Boxplot for Skewed Distributions', \emph{Computational Statistics and Data Analysis}, 52, pp. 5186-5201.

}
\author{
Marcello D'Orazio \email{mdo.statmatch@gmail.com}
}


\seealso{
\code{\link[robustbase]{adjboxStats}}, \code{\link[Hmisc]{wtd.quantile}} 
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

out <- boxB(x = x, k = 1.5, method = 'adjbox')
out$fences
out$outliers
x[out$outliers]

x[24] <- NA
x.ids <- paste0('obs',1:30)
out <- boxB(x = x, k = 1.5, method = 'adjbox', id = x.ids)
out$excluded
out$fences
out$outliers

set.seed(111)
w <- round(runif(n = 30, min=1, max=10))
out <- boxB(x = x, k = 1.5, method = 'adjbox', id = x.ids, weights = w)
out$excluded
out$fences
out$outliers

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ survey }
