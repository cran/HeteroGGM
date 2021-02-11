#' The summary of the resulting network structures.
#'
#' @name summary_network
#' @usage summary_network(opt_Mu_hat, opt_Theta_hat, data)
#' @description Summarize the characteristics of the resulting network structures.
#' @param opt_Mu_hat A p * K0_hat matrix, the optional mean vectors of K0_hat subgroups.
#' @param opt_Theta_hat n * p * K0_hat matrix, the optional precision matrices of K0_hat subgroups.
#' @param data A n * p matrix, the design matrix.
#'
#' @return A list including the overlap of edges of different subgroups, the number of edges, and the names of connected nodes to each nodes in each subgroup.
#' @export
#'
#' @import igraph
#'
#'
summary_network <- function(opt_Mu_hat, opt_Theta_hat, data){

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: summary_network
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Summarize the characteristics of the resulting network structures.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: edge_index()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ opt_Mu_hat: the optional mean vectors of K0_hat subgroups.
  ## @ opt_Theta_hat, the optional precision matrices of K0_hat subgroups.
  ## @ data: n * p matrix, the design matrix.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including:
  ## @ Theta_summary: A list including:
  ##                  @ overlap: the overlap of edges of different subgroups.
  ##                  @ Theta_nonzero_num: the number of edges of each subgroup.
  ##                  @ network_edge_summary: the names of connected nodes to each node in each subgroup.
  ## @ Mu_summary: A list including:
  ##                  @ nonzero_num: non-zero variable numbers in each subgroup.
  ##                  @ nonzero_names: non-zero variable names in each subgroup.
  ## @ variable_names: A vector, names of variables.
  ## ------------------------------------------------------------------------------------------------------------------------------------------

  variable_names <- names(data)
  if(length(variable_names) == 0){variable_names <- as.character(c(1:dim(data)[2]))}
  p <- dim(opt_Theta_hat)[1]
  K_hat <- dim(opt_Theta_hat)[3]

  if(K_hat == 1){print("warning: only one cluster.")}else{
    Mu_summary <- list()
    for (k in 1:K_hat) {
      non_k <- which(opt_Mu_hat[k,]!=0)
      if(length(non_k) > 0){
        nonzero_num <- non_k
        nonzero_names <- variable_names[nonzero_num]
        Mu_nonzero_info <- as.data.frame(cbind(nonzero_num,nonzero_names))
        Mu_summary[[k]] <- Mu_nonzero_info
      } else {
        Mu_summary[[k]] <- "ALL ZERO"
      }

    }

    Theta_nonzero_num <- list()
    for (kk in 1:K_hat) {
      Theta_nonzero_num[[kk]] <- which(opt_Theta_hat[,,kk]!=0)
    }
    overlap <- as.data.frame(matrix(0,K_hat,K_hat))
    names(overlap) <- paste("subgroup",1:K_hat)
    row.names(overlap) <- paste("subgroup",1:K_hat)
    for (k in 1:K_hat) {
      edge1 <- Theta_nonzero_num[[k]]
      for (kk in 1:K_hat) {
        edge2 <- Theta_nonzero_num[[kk]]
        if(k!=kk){
          overlap[k,kk] <- length(intersect(edge1,edge2)) - p
          overlap[kk,k] <- length(intersect(edge1,edge2)) - p
        } else{
          overlap[k,kk] <- length(intersect(edge1,edge2)) - p
          overlap[kk,k] <- length(intersect(edge1,edge2)) - p
        }
      }
    }
    network_edge_summary <- edge_index(opt_Theta_hat, data)
    Theta_summary <- list(overlap=overlap,
                          Theta_nonzero_num=apply(opt_Theta_hat, 3, function(a){sum(a!=0)}),
                          network_edge_summary=network_edge_summary)

    return(list(Theta_summary=Theta_summary, Mu_summary=Mu_summary,
                variable_names=variable_names, opt_Theta_hat=opt_Theta_hat))
  }

}

