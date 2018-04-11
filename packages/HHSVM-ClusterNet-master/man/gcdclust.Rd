\name{gcdclust}
\alias{gcdclust}

\title{gcdclust is used for calculating the solutions paths of HHSVM-CEN

}
\description{
This function is used to calculat the solutions of HHSVM-CEN 
}
\usage{
gcdclust(x, y, KK, nlambda = 30, method = "hhsvm", lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04), lambda = NULL, lambda2 = 0, pf = rep(1, nvars), pf2 = rep(1, nvars), exclude, dfmax = nvars + 1, pmax = min(dfmax * 1.2, nvars), standardize = TRUE, eps = 1e-05, maxit = 1e+06, delta = 2)
}

\arguments{
  \item{x}{
The matrix or data 
}
  \item{y}{
A response 
}
  \item{KK}{
number of cluster in data
}
  \item{nlambda}{
The number of lambda1 used in the algorithm, the defaut value is 30 
}
  \item{method}{
Our method HHSVM-CEN
}
  \item{lambda.factor}{
ifelse(nobs < nvars, 0.01, 1e-04), depends of data (relation betwin nobs and nvars)
}
  \item{lambda}{
A grid of lambda1 is calculating by our algorithm and we can use an other grid of lambda if we want.
}
  \item{lambda2}{
The defaut value is 0.
}
  \item{pf}{
 A vector if we want penalized some predector 
}
  \item{pf2}{
A vector if we want penalized some predector 
}
  \item{exclude}{
%%     ~~Describe \code{exclude} here~~
}
  \item{dfmax}{
%%     ~~Describe \code{dfmax} here~~
}
  \item{pmax}{
%%     ~~Describe \code{pmax} here~~
}
  \item{standardize}{
%%     ~~Describe \code{standardize} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The result is object of 9 items.
}
\references{
HHSVM-CEN
}
\author{

}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(test)
x=test$x
y=test$y
m=gcdclust(x, y,KK=3,lambda2=1)

}



