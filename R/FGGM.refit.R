#' Refitting of FGGM
#'
#' @author Mingyang Ren, Sanguo Zhang, Qingzhao Zhang, Shuangge Ma. Maintainer: Mingyang Ren <renmingyang17@mails.ucas.ac.cn>.
#' @references Ren, M., Zhang S., Zhang Q. and Ma S. (2020). Gaussian Graphical Model-based Heterogeneity Analysis via Penalized Fusion. Biometrics, Published Online.
#' @usage FGGM.refit(data, K, lambda1 = 0.5, lambda2 = 0.2, lambda3 = 2, a = 3, rho = 1,
#'                   eps = 5e-2, niter = 20, maxiter=10, maxiter.AMA=5,
#'                   initialization=TRUE, initialize, average=FALSE,
#'                   asymmetric=TRUE, local_appro=TRUE, penalty = "MCP", theta.fusion=TRUE)
#' @description Refitting when K0 is identified using FGGM().
#' @param data n * p matrix, the design matrix.
#' @param K Int, a selected upper bound of K_0.
#' @param lambda1 A float value, the tuning parameter controlling the sparse of the mean parameter.
#' @param lambda2 A float value, the tuning parameter controlling the sparse of the precision matrix.
#' @param lambda3 A float value, the tuning parameter controlling the number of subgroup.
#' @param a A float value, regularization parameter in MCP, the default setting is 3.
#' @param rho A float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
#' @param eps A float value, algorithm termination threshold.
#' @param niter Int, maximum number of cycles of the algorithm, the default setting is 20.
#' @param maxiter Int, maximum number of cycles of the ADMM algorithm.
#' @param maxiter.AMA Int, maximum number of cycles of the AMA algorithm.
#' @param initialization The logical variable, whether to calculate the initial value, the default setting is T, if initialization = F, the initial value uses initialize.
#' @param initialize A given initial value used if initialization = F.
#' @param average The logical variable, whether to use averaging when integrating parameters that are identified as identical subgroups, the default setting is F, which means the estimated parameters for the subgroup with the largest sample size among the subgroups identified as identical subgroups is used as the final parameter for this subgroup.
#' @param asymmetric The logical variable, symmetry of the precision matrices or not, the default setting is T.
#' @param local_appro The logical variable, whether to use local approximations when updating mean parameters, the default setting is T.
#' @param penalty The type of the penalty, which can be selected from c("MCP", "SCAD", "lasso").
#' @param theta.fusion Whether or not the fusion penalty term contains elements of the precision matrices. The default setting is T.
#'
#' @return A list including all estimated parameters and the BIC value after refitting.
#' @export
#'
#'
FGGM.refit = function(data, K, lambda1 = 0.5, lambda2 = 0.2, lambda3 = 2, a = 3, rho = 1,
                      eps = 5e-2, niter = 20, maxiter=10, maxiter.AMA=5, initialization=TRUE, initialize,
                      average=FALSE, asymmetric=TRUE, local_appro=TRUE, penalty = "MCP", theta.fusion=TRUE){

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: FGGM.refit
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Refitting when K0_hat is identified.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: FGGM()    Update_Theta()    AMA_XI()
  ##            R packages: NbClust
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ K: int, a selected upper bound of K_0.
  ## @ lambda1: a float value, the tuning parameter controlling the sparse of the mean parameter.
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the number of subgroup.
  ## @ a: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ rho: a float value, the penalty parameter in ADMM algorithm of updating precision matrix Theta, the default setting is 1.
  ## @ eps: a float value, algorithm termination threshold.
  ## @ niter: int, Maximum number of cycles of the algorithm, the default setting is 20.
  ## @ maxiter: int, Maximum number of cycles of the ADMM algorithm.
  ## @ maxiter.AMA: int, Maximum number of cycles of the AMA algorithm.
  ## @ initialization: the logical variable, whether to calculate the initial value, the default setting is T,
  ##                                     if initialization = F, the initial value uses initialize.
  ## @ initialize: A given initial value used if initialization = F.
  ## @ average: the logical variable, whether to use averaging when integrating parameters that are identified as identical subgroups,
  ## @          the default setting is F, which means the estimated parameters for the subgroup with the largest sample size among
  ## @          the subgroups identified as identical subgroups is used as the final parameter for this subgroup.
  ## @ asymmetric: the logical variable, symmetry of the precision matrices or not, the default setting is T.
  ## @ local_appro: the logical variable, whether to use local approximations when updating mean parameters, the default setting is T.
  ## @ theta.fusion: Whether or not the fusion penalty term contains elements of the precision matrices.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list PP.refit including:
  ## @ mu: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups (K0_hat is the number of identified subgroups).
  ## @ Theta: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ Xi: p * p * K0_hat array, the estimated dual variables of the precision matrices in ADMM.
  ## @ niter: int, the actual number of cycles of the algorithm.
  ## @ diff_Xi: a float value, the L2-norm difference between subgroups.
  ## @ prob0: K0_hat * 1 vector, the estimated mixture probabilities of subgroups.
  ## @ L.mat0: n * K0_hat matrix, the estimated probability that each sample belongs to each subgroup.
  ## @ Theta0: p * p * K array, the estimated original precision matrices of given K subgroups.
  ## @ Xi0: p * p * K array, the estimated original dual variables of the precision matrices of given K subgroups in ADMM.
  ## @ mu0: K * p matrix, the estimated original mean vectors of K subgroups.
  ## @ group: a list, the sequence number partition of original K subgroups compressed to K0_hat.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## @ fit.error: a float value, the value of the loss function (without penalty function) corresponding the choice of given tuning parameters.
  ## @ df: a float value, the penalty value for non-zero parameters corresponding the choice of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------

  PP = FGGM(data, K, lambda1, lambda2, lambda3, eps=eps, maxiter=maxiter, maxiter.AMA=maxiter.AMA,
            initialization=initialization, initialize=initialize, average=average,
            asymmetric=asymmetric, local_appro=local_appro, penalty = penalty, theta.fusion=theta.fusion)
  K_hat = length(PP$group)
  if(K_hat == 1){
    return(PP)
  } else {
    PP.refit = FGGM(data, K_hat, lambda1, lambda2, 0, eps=eps, maxiter=maxiter, maxiter.AMA=maxiter.AMA,
                    initialization=TRUE, initialize=initialize, average=average,
                    asymmetric=asymmetric, local_appro=local_appro, penalty = penalty, theta.fusion=theta.fusion)
  }
  return(PP.refit)
}
