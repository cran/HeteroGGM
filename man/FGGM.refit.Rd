% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FGGM.refit.R
\name{FGGM.refit}
\alias{FGGM.refit}
\title{Refitting of FGGM}
\usage{
FGGM.refit(data, K, lambda1 = 0.5, lambda2 = 0.2, lambda3 = 2, a = 3, rho = 1,
                  eps = 5e-2, niter = 20, maxiter=10, maxiter.AMA=5,
                  initialization=TRUE, initialize, average=FALSE,
                  asymmetric=TRUE, local_appro=TRUE, penalty = "MCP", theta.fusion=TRUE)
}
\arguments{
\item{data}{n * p matrix, the design matrix.}

\item{K}{Int, a selected upper bound of K_0.}

\item{lambda1}{A float value, the tuning parameter controlling the sparse of the mean parameter.}

\item{lambda2}{A float value, the tuning parameter controlling the sparse of the precision matrix.}

\item{lambda3}{A float value, the tuning parameter controlling the number of subgroup.}

\item{a}{A float value, regularization parameter in MCP, the default setting is 3.}

\item{rho}{A float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.}

\item{eps}{A float value, algorithm termination threshold.}

\item{niter}{Int, maximum number of cycles of the algorithm, the default setting is 20.}

\item{maxiter}{Int, maximum number of cycles of the ADMM algorithm.}

\item{maxiter.AMA}{Int, maximum number of cycles of the AMA algorithm.}

\item{initialization}{The logical variable, whether to calculate the initial value, the default setting is T, if initialization = F, the initial value uses initialize.}

\item{initialize}{A given initial value used if initialization = F.}

\item{average}{The logical variable, whether to use averaging when integrating parameters that are identified as identical subgroups, the default setting is F, which means the estimated parameters for the subgroup with the largest sample size among the subgroups identified as identical subgroups is used as the final parameter for this subgroup.}

\item{asymmetric}{The logical variable, symmetry of the precision matrices or not, the default setting is T.}

\item{local_appro}{The logical variable, whether to use local approximations when updating mean parameters, the default setting is T.}

\item{penalty}{The type of the penalty, which can be selected from c("MCP", "SCAD", "lasso").}

\item{theta.fusion}{Whether or not the fusion penalty term contains elements of the precision matrices. The default setting is T.}
}
\value{
A list including all estimated parameters and the BIC value after refitting.
}
\description{
Refitting when K0 is identified using FGGM().
}
\references{
Ren, M., Zhang S., Zhang Q. and Ma S. (2020). Gaussian Graphical Model-based Heterogeneity Analysis via Penalized Fusion. Biometrics, Published Online.
}
\author{
Mingyang Ren, Sanguo Zhang, Qingzhao Zhang, Shuangge Ma. Maintainer: Mingyang Ren <renmingyang17@mails.ucas.ac.cn>.
}
