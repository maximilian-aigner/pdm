\name{coef.gcdclust}
\alias{coef.gcdclust}
\title{get coefficients or make coefficient predictions from a "gcdclust" object.}
\description{
This function gets coefficients or makes coefficient predictions from a gcdclust model,
using the stored \code{"gcdclust.fit"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{coef}{gcdclust}(object,s=c("lambda.1se","lambda.min"),...)
}
\arguments{
	\item{object}{fitted \code{\link{gcdclust}} object.}
	\item{s}{value(s) of the penalty parameter \code{lambda} at which
	predictions are required. Default is the value \code{s="lambda.1se"} stored
	on the \code{object}, it is the largest value of \code{lambda} such that error is
	within 1 standard error of the minimum. Alternatively \code{s="lambda.min"} can be
	used, it is the optimal value of \code{lambda} that gives minimum
	cross validation error \code{cvm}. If \code{s} is numeric, it is taken as the value(s) of
	\code{lambda} to be used.}
	\item{\dots}{not used. Other arguments to predict. } 
}
\details{This function makes it easier to use the results of
a gcdclust model to get coefficients or make coefficient predictions.}
\value{The object returned depends the \dots argument which is passed on
to the \code{\link{predict}} method for \code{\link{gcdclust}} objects.}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
BugReport: \url{http://code.google.com/p/gcdclust/}\cr

Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\seealso{\code{\link{cv.gcdclust}}, and \code{\link{predict.cv.gcdclust}} methods.}
\examples{

data(test)
x=test$x
y=test$y
m=gcdclust(x, y,KK=3,lambda2=1)
coef=coef(m,type="coef")

#NextMethod("coef")
}

