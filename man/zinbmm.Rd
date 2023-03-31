\name{zinbmm}
\alias{zinbmm}
\title{Zero-Inflated Negative Binomial mixture model for clustering and gene selection}
\usage{
zinbmm(mat,Bmat,lib=NULL,K,phi_init=NULL,phi_global=T,pi_ind=F,
       tune_set=NULL,ntune=10,c='BIC',ncores=1,maxit=500,tol=1e-05,vrbs=F)
}
\description{
The main function of the ZINBMM package. Applies simultaneous clustering and cluster-specific gene selection for scRNA-seq raw data while accounting for the potential dropouts and batch effects.
}
\arguments{
\item{mat}{matrix of expression count data, with each row representing one sample/cell.}
\item{Bmat}{defines the batch information of each sample/cell.}
\item{lib}{defines the library size of each sample/cell.}
\item{K}{specifies thenumber of the clusters.}
\item{phi_init}{a vector of the initialized values of the dispersion parameters, default to \code{NULL}.}
\item{phi_global}{a boolean variable which defines whether the dispersion parameters are fixed across genes, default to \code{TRUE}.}
\item{pi_ind}{a boolean variable which defines whether the zero inflation parameters are independent of the clusters, default to \code{FALSE}.}
\item{tune_set}{optional user-supplied tuning sequence; default is \code{NULL}, and ZINBMM chooses its own sequence.}
\item{ntune}{specifies the number of the tuning parameters with default value \code{10}.}
\item{c}{criterion to determine the optimal tuning parameter. The default is \code{c=='BIC'}, which uses the Bayesian information criterion. Also allows for \code{c=='AIC'} and \code{c=='mBIC'} for Akaike information criterion and modified BIC.}
\item{ncores}{defines the number of cores to be used on the user's machine, default to \code{1}.}
\item{maxit}{defines the maximum number of iterations, default to \code{500}.}
\item{tol}{specifies the iteration convergency tolerance, default to \code{1e-05}.}
\item{vrbs}{a boolean variable which defines whether to print the iteration details, default to \code{FALSE}.}}

\value{
The \code{zinbmm} returns a list containing the following components:
\item{clust}{the consecutive integers specifying the clustering index of each sample/cell.}
\item{pclust}{the mixing proportions of \code{K} clusters.}
\item{beta_info}{the indicators of the selected informative genes.}
\item{batch_effect}{the estimates for batch effect parameters.}
\item{zero_infl}{the estimates for zero-inflated parameters.}
}


\examples{
# This is an example of zinbmm with simulated data
data(example)
mat=example$mat
Bmat=example$batch_mat
# Apply zinbmm
results=zinbmm(mat,Bmat,K=3,ntune=8,tol=1e-03,ncores=2)
## clustering
results$clust
## gene selection
results$beta_info

}
