\name{plot4ratios}
\alias{plot4ratios}
\title{
Draws a scatter-plot that summarizes the findings of the Hidiroglou-Berthelot outliers' or the ratioSize detection method.
}

\description{
The function gets the output of the function of \code{HBmethod} or \code{ratioSize} when they are ran with the argument \code{return.dataframe = TRUE}) to draw a scatter-plot of ratios vs. the corresponding importance measures.  
}

\usage{
plot4ratios(out)
}
\arguments{
  \item{out}{
Is the output of  \code{\link[univOutl]{HBmethod}} or \code{\link[univOutl]{ratioSize}} when they are launched with the argument \code{return.dataframe = TRUE}.
}
}

\details{
This function draws a scatter-plot. With the output of \code{HBmethod} the ratios (=yt2/yt1) are on Y axis while their importance measure ( max(yt1, yt2)^U) are represented on the X axis. With the output of \code{ratioSize} on the Y axis the centered ratios are reported. In addition the acceptance bounds are drawn (blue lines); the dots (in red color) outside the bounds are the outliers. This is considered a useful diagnostic plot to understand how the procedure identifies the outliers.  
}

\value{
A scatter-plot is drawn and, in addition, the output includes a data.frame with the data used to derive the plot.
}
 
\references{
Hidiroglou, M.A. and Berthelot, J.-M. (1986) `Statistical editing and Imputation for Periodic Business Surveys'. \emph{Survey Methodology}, Vol 12, pp. 73-83.

Hidiroglou, M.A. and Emond, N. (2018) `Modifying the Hidiroglou-Berthelot (HB) method'. \emph{Unpublished note}, Business Survey Methods Division, Statistics Canada, May 18 2018.
}

\author{
Marcello D'Orazio \email{mdo.statmatch@gmail.com}
}

\seealso{
\code{\link[univOutl]{HBmethod}}, \code{\link[univOutl]{ratioSize}}
}


\examples{
# generate some data
set.seed(222)
x0 <- rnorm(30, 50, 5)
set.seed(333)
rr <- runif(30, 0.9, 1.2)
rr[10] <- 2
x1 <- x0 * rr

# run HBmethod with argument return.dataframe = TRUE
out <- HBmethod(yt1 = x0, yt2 = x1,  
                return.dataframe = TRUE)
# draw the scatterplot
plot4ratios(out)

}
\keyword{ univar }
\keyword{ survey }
