\name{plot.gcdclust}
\alias{plot.gcdclust}
\title{Plot coefficients from a "gcdclust" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{gcdclust}} object. This function is modified based on the \code{plot} function from the \code{glmnet} package.
}
\usage{
\method{plot}{gcdclust}(x, xvar = c("norm", "lambda"), color = FALSE, label = FALSE, ...)
}
\arguments{
		\item{x}{fitted \code{\link{gcdclust}} model}
		\item{xvar}{what is on the X-axis. \code{"norm"} plots against the
		L1-norm of the coefficients, \code{"lambda"} against the log-lambda
		sequence.}
		\item{color}{if \code{TRUE}, plot the curves with rainbow colors. \code{FALSE} is gray colors. Default is \code{FALSE}}
		\item{label}{if \code{TRUE}, label the curves with variable sequence numbers. Default is \code{FALSE}}
		\item{\dots}{other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced.
}

\references{
HHSVM-CEN
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr



Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\examples{

data(test)
x=test$x
y=test$y


m=gcdclust(x,y,KK=3,lambda2=5)
plot(m)
} 


