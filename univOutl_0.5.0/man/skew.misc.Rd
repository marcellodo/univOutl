\name{skew.misc}
\alias{skew.misc}
\title{
Calculates some skewness measures.
}

\description{
The function calculates some skewness measures for the input vector data.  
}

\usage{
skew.misc(x, weights=NULL)
}
\arguments{
  \item{x}{
Input vector containing data for which skewness will be calculated.
}
  \item{weights}{
  Optional vector with eventual non-negative weights associated to the units in \code{x}
 }
}

\details{
This function calculates Pearson's skewness coefficient, the MedCouple measure of skewness and the non-parametric Bowley's measure of symmetry. The Bowley's skewness measure uses quartiles:

\deqn{b = \frac{(Q_3 - Q_2) - (Q_2 - Q_1) }{Q_3 - Q_1}  }{b = ((Q3 - Q2) - (Q2 - Q1))/(Q3 - Q1)}

It ranges between -1 and +1, where positive (negative) values denote right (left) skewness. A value equal to 0 indicates symmetry.
A crude measure of skewness can be obtained with a monotonic increasing function of b:

\deqn{ g=\frac{1+b}{1-b} }{g = (1+b)/(1-b)}

It ranges from 0 to Inf, g=1 indicates symmetry.

A  measure of skewness similar to the Bowley's one is achieved by replacing Q3 and Q1 with respectively P90 and P10 percentiles:

\deqn{b_P = \frac{(P_{90} - P_{50}) - (P_{50} - P_{10}) }{P_{90} - P_{10}}  }{b_P = ((P90 - P50) - (P50 - p10))/(P90 - P10)}

Similarly
\deqn{ g_P=\frac{1+b_P}{1-b_P} }{g_P = (1+b_P)/(1-b_P)}

For major details see Kotz at al. (2006, vol. 12, pp. 7771-7772).

The medCouple measure of skewness, M, ranges from -1 to +1 and is equal to 0 in case of symmetry, while \eqn{M > 0}{M>0} indicates positive skewness. For major details see \code{\link[robustbase]{mc}}.

Note that eventual weights, passed through the argument \code{weights}, are used ONLY in the calculation of the Bowley's type measures. 
  
}

\value{
A vector with the estimated measures of skewness.
}
 
\references{
Kotz S. et al. (2006) \emph{Encyclopedia of Statistical Sciences, Volume 12}. John Wiley and Sons.
}

\author{
Marcello D'Orazio \email{mdo.statmatch@gmail.com}
}

\seealso{
\code{\link[robustbase]{mc}}
}


\examples{

set.seed(112233)
y <- rnorm(n = 30, mean = 50, sd = 10)
y[20] <- 100

skew.misc(x = y, weights=NULL)

# use weights
ww <- runif(n = 30, min = 1, max = 10)
skew.misc(x = y, weights=ww)

}
\keyword{ univar }
\keyword{ survey }
