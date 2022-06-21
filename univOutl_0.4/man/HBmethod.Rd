\name{HBmethod}
\alias{HBmethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hidiroglou-Berthelot procedure for detecting outliers with periodic data
}
\description{
This function implements the method proposed by Hidiroglou and Berthelot (1986) to identify outliers in periodic data, i.e. when the same variable is measured at two time points.
}
\usage{
HBmethod(yt1, yt2, U=0.5, A=0.05, C=4, pct=0.25,
         id=NULL, std.score=FALSE, return.dataframe=FALSE, adjboxE=FALSE)

}

\arguments{

  \item{yt1}{
Numeric vector providing the values observed at time t1.
}
  \item{yt2}{
Numeric vector providing the values observed at time t2 (t2 > t1).
}
  \item{U}{
Numeric, parameter needed to determine the `importance' of a ratio. The value should lie in [0, 1] interval; commonly used values are 0.3, 0.4, or 0.5 (default) (see Details for further information).
}
  \item{A}{
Numeric, parameter needed when computing the scale measure used to derive the bounds. Hidiroglou and Berthelot (1986) suggest setting \eqn{A = 0.05}{A = 0.05} (default) (see Details for further information).
}
  \item{C}{
Numeric, parameter determining the extension of the interval; greater values will provide larger intervals, i.e. fewer expected outliers. Values commonly used are 4 (default) or 7, but also values close or grater than 40 can be used in some particular cases. Note that two C values can be provided instead of one, the first one will be used to determine the left tail bound, while the second determines the right tail bound; this setting can help in improving outlier detection in skewed distributions (see Details for further information).
}
  \item{pct}{
Numeric, the percentage point of the scores that will be used to calculate the lower and upper bounds. By default, \code{pct = 0.25}, i.e. quartiles Q1 and Q3 are considered. In some cases, as suggested by Hidiroglou and Emond (2018), using \code{pct = 0.10}, i.e. percentiles P10 and P90, may be a better choice. Se Details for further information.   
}
  \item{id}{
Optional numeric or character vector, with identifiers of units. If \code{id=NULL} units identifiers will be set equal to their position.
}

  \item{std.score}{
Logical, if \code{TRUE} the output will include a standardized score variable (see Details, for further information)
}
  \item{return.dataframe}{
Logical, if \code{TRUE} the output will save all the relevant information for outlier detection in a dataframe with the following columns: 'id' (units' identifiers), `yt1', `yt2', `ratio' (= yt1/yt2), `sizeU' (=max(yt1, yt2)^U),`Escore' (the E scores, see Details), `std.Escore' (the standardized E scores when \code{std.score=TRUE}, see Details) and finally `outliers', where value 1 indicates observations detected as an outlier, 0 otherwise.  
}
 \item{adjboxE}{
Logical (default \code{FALSE}), if \code{TRUE} an additional search of outliers will be done on the E-scores using the boxplot adjusted for skewness as implemented in the function \code{\link[univOutl]{boxB}} when run with with argument \code{method = "adjbox"}.
}

}

\details{
The method proposed by Hidiroglou and Berthelot (1986) to identify outliers in periodic data consists in deriving a score variable based on the ratios \eqn{r_i = y_{i,t2}/y_{i,t1}}{r_i = y_{i,t2}/y_{i,t1}} (\code{yt2/yt1}) with \eqn{i=1,2,\ldots, n}{i=1,2,...,n} being \eqn{n}{n} the number of observations after discarding NAs and 0s in both \code{yt1} and \code{yt2}.

At first the ratios are centered around their median \eqn{r_M}{r_M}:

\deqn{s_i = 1 - r_M/r_i \qquad \mbox{if} \qquad 0 < r_i < r_M}{s_i = 1 - r_M/r_i  if  0 < r_i < r_M}

\deqn{s_i = r_i/r_M - 1  \qquad \mbox{if} \qquad r_i \geq r_M}{s_i = r_i/r_M - 1  if  r_i >= r_M}

Then, in order to account for the magnitude of data, the following score is derived:

\deqn{E_i = s_i \times [ max(y_{i,t1},y_{i,t2} ) ]^U}{E_i = s_i [max(y_{i,t1},  y_{i,t2})]^U} 

Finally, the interval is calculated as:

\deqn{(E_M-C \times d_{Q1}, E_M+C\times d_{Q3} )}{(E_M - C*d_{Q1}, E_M + C*d_{Q3})}

where 

\eqn{d_{Q1} = max( E_M - E_{Q1}, |A \times E_M| ) }{max( E_M - E_{Q1}, |A*E_M| )} and \eqn{d_{Q3} = max( E_{Q3} - E_M, |A \times E_M| ) }{max( E_{Q3} - E_M, |A*E_M| )}

being \eqn{E_{Q1}}{E_Q1}, \eqn{E_M}{E_M} and \eqn{E_{Q3}}{E_Q3} the quartiles of the E scores (when \code{pct = 0.25}, default)).
 Recently Hidiroglou and Emond (2018) suggest using percentiles P10 and P90 of the E scores in replacement of respectively Q1 and Q3  to avoid the drawback of many units identified as outliers; this is likely to occur when a large proportion of units (>1/4) has the same ratio. P10 and P90 are achieved by setting \code{pct = 0.10} when running the function. 
 
In practice, all the units with an E score outside the interval are considered as outliers. Notice that when two \code{C} values are provided, then the first is used to derive the left bound while the second determines the right bound.

When \code{std.score=TRUE} a standardized score is derived in the following manner:

\deqn{z_{E,i} = g \times \frac{E_i - E_M}{d_{Q1}} \qquad \mbox{if} \qquad  E_i < E_M}{z_E,i = g*(E_i - E_M)/d_{Q1}  if  E_i < E_M}

\deqn{z_{E,i} = g \times \frac{E_i - E_M}{d_{Q3}} \qquad \mbox{if} \qquad  E_i \geq E_M}{z_E,i = g*(E_i - E_M)/d_{Q3}  if  E_i >= E_M}

The constant g is set equal to \code{qnorm(1-pct)} and  makes \eqn{d_{Q1}}{d_{Q1}} and  \eqn{d_{Q3}}{d_{Q3}} approximately unbiased estimators when the E scores follow the normal distribution.

When \code{adjboxE = TRUE} outliers on the E scores will all be searched using the boxplot adjusted for skewness as implemented in the function \code{\link[univOutl]{boxB}} when run with with argument \code{method = "adjbox"}.
}

\value{
A list whose components depend on the \code{return.dataframe} argument. When \code{return.dataframe=FALSE} just the following components are provided:

\item{median.r}{the median of the ratios}
\item{quartiles.E}{Quartiles of the E score}
\item{bounds.E}{Bounds of the interval of the E score, values outside are considered outliers.}
\item{excluded}{
	The identifiers or positions (when \code{id=NULL}) of units in both \code{yt1} and \code{yt2} that are excluded by the outliers detection, i.e. NAs and 0s.
}
\item{outliers}{
The identifiers or positions (when \code{id=NULL}) of units in \code{yt1} or \code{yt2} identified as outliers.
}
\item{outliersBB}{
The identifiers or positions (when \code{id=NULL}) of units in \code{yt1} or \code{yt2} identified as outliers by applying the boxplot adjusted for skewness to the E scores. This component appears in the output only when \code{adjboxE = TRUE}.
}

When \code{return.dataframe=TRUE}, the first three components remain the same with, in addition, two dataframes:

\item{excluded}{
	A dataframe with the subset of observations excluded. The data frame has the following columns: 'id' (units' identifiers), 'yt1' columns 'yt2'}

\item{data}{
	A dataframe with the the not excluded observations and the following columns: `id' (units' identifiers), `yt1', `yt2', `ratio' (= yt1/yt2), `sizeU' (=max(yt1, yt2)^U),`Escore' (the E scores, see Details), `std.Escore' (the standardized E scores when \code{std.score=TRUE}, see Details) and `outliers', where value 1 indicates observations detected as an outlier, 0 otherwise. in addition the column `outliersBB' will also be included when \code{adjboxE = TRUE}. 
}

}

\references{
Hidiroglou, M.A. and Berthelot, J.-M. (1986) `Statistical editing and Imputation for Periodic Business Surveys'. \emph{Survey Methodology}, Vol 12, pp. 73-83.

Hidiroglou, M.A. and Emond, N. (2018) `Modifying the Hidiroglou-Berthelot (HB) method'. \emph{Unpublished note}, Business Survey Methods Division, Statistics Canada, May 18 2018.

}
\author{
Marcello D'Orazio \email{mdo.statmatch@gmail.com}
}

\seealso{
\code{\link[univOutl]{plot4ratios}}, \code{\link[univOutl]{ratioSize}}
}

\examples{

set.seed(222)
x0 <- rnorm(30, 50, 5)
x0[1] <- NA
set.seed(333)
rr <- runif(30, 0.9, 1.2)
rr[10] <- 2
x1 <- x0 * rr
x1[20] <- 0

out <- HBmethod(yt1 = x0, yt2 = x1)
out$excluded
out$median.r
out$bounds.E
out$outliers
cbind(x0[out$outliers], x1[out$outliers])

out <- HBmethod(yt1 = x0, yt2 = x1,  
                return.dataframe = TRUE)
out$excluded
head(out$data)

}
\keyword{ univar }
\keyword{ survey }

