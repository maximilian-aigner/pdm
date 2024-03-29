% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastphase.R
\name{SNPknock.fp.loadFit}
\alias{SNPknock.fp.loadFit}
\title{Load the parameter estimates obtained by fastPhase and assembles the HMM model for the genotype data.}
\usage{
SNPknock.fp.loadFit(r_file, theta_file, alpha_file, x)
}
\arguments{
\item{r_file}{a string with the path of the "_rhat.txt" file produced by fastPhase.}

\item{theta_file}{a string with the path of the "_thetahat.txt" file produced by fastPhase.}

\item{alpha_file}{a string with the path of the "_alphahat.txt" file produced by fastPhase.}

\item{x}{a numpy array of length p, where p is the number of SNPs, containing the genotype sequence of the first individual in the dataset.}
}
\value{
A structure describing the HMM fitted by fastPhase.
}
\description{
This function loads the parameter estimates obtained by fastPhase (see \link{SNPknock.fp.runFastPhase})
and assembles the HMM model for the genotype data, in the format required by the knockoff generation function
\link{SNPknock.knockoffHMM}.
}
\details{
In the description of the parameter "x", the first individual is intended 
in the same order as provided to fastPhase. This is needed in order to correctly 
intepret the emission parameters estimated by fastPhase.

This function returns a structure with three fields: 
\itemize{
  \item{"pInit": an array of length K, containing the marginal distribution of the hidden states for the first SNP.}
  \item{"Q": an array of size (p-1,K,K), containing a list of p-1 transition matrices between the K latent states of the HMM.}
  \item{"pEmit": an array of size (p,K,3), containing the emission probabilities of the hidden states for each of the p SNPs.}
 }
}
\examples{
# Load an example data matrix X from the package installation directory.
X_file = system.file("extdata", "X.RData", package = "SNPknock")
load(X_file)

# Specify the location of the fastPhase output files containing the parameter estimates.
# Example files can be found in the package installation directory.
r_file = system.file("extdata", "X_rhat.txt", package = "SNPknock")
theta_file = system.file("extdata", "X_thetahat.txt", package = "SNPknock")
alpha_file = system.file("extdata", "X_alphahat.txt", package = "SNPknock")

# Read the parameter files and build the HMM
hmm = SNPknock.fp.loadFit(r_file, theta_file, alpha_file, X[1,])

}
\references{
Scheet and Stephens,  A fast and flexible statistical model for large-scale population genotype data,
  Am J Hum Genet (2006).
  \href{http://www.sciencedirect.com/science/article/pii/S000292970763701X}{http://www.sciencedirect.com/science/article/pii/S000292970763701X}
}
\seealso{
Other fastphase: \code{\link{SNPknock.fp.runFastPhase}},
  \code{\link{SNPknock.fp.writeX}}
}
