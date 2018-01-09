\name{LocScaleB}
\alias{LocScaleB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate outlier detection with bounds based on robust location and scale estimates
}
\description{
This function identifies outliers in the tails of a distribution by detecting the observations outside the bounds built using a robust estimate of both location and scale parameters.  
}
\usage{
LocScaleB(x, k=3, method='MAD',  weights=NULL, id=NULL, 
          exclude=NA, logt=FALSE, return.dataframe=FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
		Numeric vector that will be searched for outliers.
}
  \item{k}{
		Nonnegative constant that determines the extension of bounds. Commonly used values are 2, 2.5 and 3 (default).
}
  \item{method}{
		character identifying how to estimate the scale of the distribution. Available choices are:
		
   	   \code{method='IQR'} for using the Inter-Quartile Range, i.e. Q75-Q25; 
 
  	   \code{method='IDR'} for using the Inter-Decile Range; i.e. Q90-Q10 
		 		   	   
		\code{method='MAD'} for using the Median Absolute Deviation;

		\code{method='Gini'} robust scale estimate based on Gini's Mean Difference (see  \code{\link[Hmisc]{GiniMd}}); 
		
		\code{method='ScaleTau2'} robust tau-estimate of univariate scale, as proposed by Maronna and Zamar (2002) (see also\code{\link[robustbase]{scaleTau2}}); 
		
		\code{method='Qn'} for using the Qn estimator proposed by Rousseeuw and Croux (1993) (see also \code{\link[robustbase]{Qn}}); 
		
		\code{method='Sn'} for using the Sn estimator proposed by Rousseeuw and Croux (1993) (see also \code{\link[robustbase]{Sn}}). 
		
		When \code{method='dQ'} the estimated scale for the left tail is (Q2-Q1)/0.6745, while for the right tail it is considered (Q3-Q2)/0.6745; this double estimate should be able to account for slight skewness.
		
		 Finally, when \code{method='AdjOut'}, bounds are based on the adjusted outlyingness method as proposed by Hubert and Van der Veeken (2008).
}

  \item{weights}{
		Optional numeric vector that provides weights associated to observations. Only nonnegative weights  are allowed. Note that weights can only be used when \code{method='MAD'}, \code{method='IQR'} or \code{method='dQ'}.
}
  \item{id}{
Optional numeric or character vector, with identifiers of units in \code{x}. If \code{id=NULL} (deafult) units' identifiers will be set equal to their position in \code{x}.
}

  \item{exclude}{
	Values of \code{x} that will be excluded by the analysis. By default missing values (\code{exclude = NA})
}
  \item{logt}{
Logical, if \code{TRUE}, before searching outliers the \code{x} variable is log-transformed (\code{log(x+1)} is considered). Note that in this case that summary output (bounds, etc.) will refer to log-transformed variable.
}
  \item{return.dataframe}{
Logical, if \code{TRUE} the output will save all the relevant information for outlier detection in a dataframe with the following columns: `id' (units' identifiers), `x', `log.x' (only if \code{logt=TRUE}), `weight' (only when argument \code{weights} is provided), `score' (the standardized scores, see Details) and, finally, `outliers', where value 1 indicates observations detected as an outlier, 0 otherwise.  
}



}
\details{
The intervals are derived by considering the median as a robust location estimate while different robust scale estimators are considered:

\deqn{[Q_2 - k \times \tilde{s}_L; \quad Q_2 + k \times \tilde{s}_R]}{ [Q2 - k*s_L; Q2  + k*s_R] }

\eqn{Q_2}{Q2} is the median, while \eqn{\tilde{s}_L}{ s_L} and \eqn{ \tilde{s}_R}{S_R} are robust scale estimates. 
With most of the methods \eqn{\tilde{s}_L=\tilde{s}_R}{ s_L=s_L} with exception of \code{method='dQ'} where:

\deqn{\tilde{s}_L = (Q_2 - Q_1)/0.6745 \qquad \mbox{and} \qquad \tilde{s}_R = (Q_3 - Q_2)/0.6745}{s_L = (Q2 - Q1)/0.6745  and s_R = (Q3 - Q2)/0.6745}

With \code{method='AdjOut'} the following estimates are considered: 

\deqn{\tilde{s}_L = (Q_2 - f_L) \qquad \mbox{and} \qquad \tilde{s}_R = (f_R - Q_2)}{s_L = (Q2 - fL) ands_R = (fR - Q2)}

being \eqn{f_R}{fR} and \eqn{f_L}{fL} derived starting from the fences of the adjusted boxplot (Hubert and Vandervieren, 2008; see \code{\link[robustbase]{adjboxStats}}). 

When weights are available (passed via the argument \code{weights}) then they are used in the computation of the quartiles. In particular, the quartiles are derived using the function  \code{\link[Hmisc]{wtd.quantile}} in the package \pkg{Hmisc}. Note that their use is allowed just with \code{method='IQR'}, \code{method='dQ'} or \code{method='AdjOut'}.

The `score' variable reported in the the \code{data} dataframe when \code{return.dataframe=TRUE} is the standardized score derived as (x - Median)/scale.  
}

\value{
A list whose components depend on the \code{return.dataframe} argument. When \code{return.dataframe = FALSE} just the following components are provided:

\item{pars}{Vector with estimated median and scale parameters}
\item{bounds}{The bounds of the interval, values outside the interval are considered outliers.}
\item{excluded}{
	The position or identifiers of \code{x} values excluded by outlier detection, according to the argument \code{exclude}
}
\item{outliers}{
	The position or identifiers of \code{x} values detected as outliers (outside bounds).
}

When \code{return.dataframe=TRUE} the latter two components are substituted with two dataframes:

\item{excluded}{
	A dataframe with the subset of observations excluded.
}

\item{data}{
	A dataframe with the the not excluded observations and the following columns: `id' (units' identifiers), `x', `log.x' (only if \code{logt=TRUE}), `weight' (only when argument \code{weights} is provided), `score' (the standardized scores, see Details) and, finally, `outliers', where value 1 indicates observations detected as an outlier and 0 otherwise.
}
}
\references{

Hubert, M. and Van der Veeken, S. (2008) `Outlier Detection for Skewed Data'. \emph{Journal of Chemometrics}, 22, pp. 235-246.

Maronna, R.A. and Zamar, R.H. (2002) `Robust estimates of location and dispersion of high-dimensional datasets' \emph{Technometrics}, 44, pp. 307-317.

Rousseeuw, P.J. and Croux, C. (1993) `Alternatives to the Median Absolute Deviation', \emph{Journal of the American Statistical Association} 88, pp. 1273-1283.

Vanderviere, E. and Huber, M. (2008) `An Adjusted Boxplot for Skewed Distributions', \emph{Computational Statistics & Data Analysis}, 52, pp. 5186-5201

}
\author{
Marcello D'Orazio \email{mdo.statmatch@gmail.com}
}


\seealso{
\code{\link[stats]{mad}}, \code{\link[robustbase]{scaleTau2}}, \code{\link[robustbase]{Qn}}, \code{\link[robustbase]{Sn}}, \code{\link[Hmisc]{GiniMd}}
}
\examples{

set.seed(333)
x <- rnorm(30, 50, 1)
x[10] <- 1
x[20] <- 100

out <- LocScaleB(x = x,  k = 3, method='MAD')
out$pars
out$bounds
out$outliers
x[out$outliers]

out <- LocScaleB(x = x,  k = 3, method='MAD',
                 return.dataframe = TRUE)
head(out$data)

out <- LocScaleB(x = x, k = 3, method='AdjOut')
out$outliers


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ robust }
\keyword{ survey }