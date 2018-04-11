\name{cv.gcdclust}
\alias{cv.gcdclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Do a cross validation 

}
\description{
This function is used to choce the best lambda.min that gives the minimum error of classification
}
\usage{
cv.gcdclust(x, y, KK, lambda = NULL, pred.loss = c("misclass", "loss"), nfolds = 5, foldid, delta = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{fitted \code{\link{cv.gcdclust}} object.}
	\item{s}{value(s) of the penalty parameter \code{lambda} at which
	predictions are required. Default is the value \code{s="lambda.1se"} stored
	on the CV \code{object}, it is the largest value of \code{lambda} such that error is
	within 1 standard error of the minimum. Alternatively \code{s="lambda.min"} can be
	used, it is the optimal value of \code{lambda} that gives minimum
	cross validation error \code{cvm}. If \code{s} is numeric, it is taken as the value(s) of
	\code{lambda} to be used.}
	 
}

\details{This function makes 
cross-validation to get lambda.min.}
\value{Result of cross validation}

\references{
HHSVM-CEN
}
\author{

}
\note{

}



\examples{
data(test)
x=test$x
y=test$y
cv=cv.gcdclust(x,y,KK=3,lambda2=1)


}

