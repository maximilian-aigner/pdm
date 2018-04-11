\name{plot.cv.gcdclust}
\alias{plot.cv.gcdclust}
\title{plot the cross-validation curve produced by cv.gcdclust}
\description{
Plots the cross-validation curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used. This function is modified based on the \code{plot.cv} function from the \code{glmnet} package.}
\usage{
\method{plot}{cv.gcdclust}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{fitted \code{\link{cv.gcdclust}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or
	its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}
\details{A plot is produced.}

\references{
HHSVM - CEN
An efficient algorithm of calculating the solutions paths of hhsvm
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.gcdclust}}.}
\examples{
data(test)
x=test$x
y=test$y
cv=cv.gcdclust(x,y,KK=3,lambda2=1)
plot(cv)

}


