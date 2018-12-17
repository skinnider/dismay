#' Calculate the estimator of Kendall's tau for zero-inflated count data
#' proposed by Pimentel (Stat Prob Lett (2015) 96:61-67) between each pair of 
#' columns in a matrix.  
#' 
#' @param a matrix of data 
#' @return the estimated correlation coefficients
#' @export
kendall_zi <- function(mat) {
  # notation from Pimentel et al. Stat Prob Lett (2015) 96:61-67 
  n <- nrow(mat)
  # calculate Kendall's tau for non-zero pairs only 
  message("calculating Kendall's tau for paired non-zero observations...")
  t11 <- cor.fk.nz(mat)
  # calculate the remaining values needed for the estimator
  message("calculating estimator of Kendall's tau for zero-inflated data...")
  nz <- mat != 0
  p11 <- crossprod(nz) / n
  p00 <- crossprod(!nz) / n
  p01 <- (t(!nz) %*% nz) / n
  p10 <- (t(nz) %*% !nz) / n
  # conditional distributions
  message("calculating conditional distribution 1 of 2...")
  p1 <- p1.mat(mat)
  message("calculating conditional distribution 2 of 2...")
  p2 <- p2.mat(mat)
  # estimator
  tHat <- p11^2 * t11 + 2 * ((p00 * p11) - (p01 * p10)) + 
    2 * p11 * (p10 * (1 - 2 * p1) + p01 * (1 - 2 * p2))
  return(tHat)
}

p1.mat <- function(mat) {
  p1 <- function(x, y) {
    y <- y[x != 0]
    x <- x[x != 0]
    x10 <- x
    x10[y != 0] <- 0
    x11 <- x
    x11[y == 0] <- 0
    mean(x10 > x11)
  }
  # get all permutations of columns
  permutations <- gtools::permutations(n = ncol(mat), r = 2, 
                                       v = seq_len(ncol(mat)))
  # create a matrix to store the results
  result <- matrix(0, nrow = ncol(mat), ncol = ncol(mat))
  # apply the function to each combination
  result[permutations] <- pbapply::pbapply(permutations, 1, function(r) p1(
    mat[, r[1]], mat[, r[2]]))
  return(result)
}

p2.mat <- function(mat) {
  p2 <- function(x, y) {
    x <- x[y != 0]
    y <- y[y != 0]
    y01 <- y
    y01[x != 0] <- 0
    y11 <- y
    y11[x == 0] <- 0
    mean(y01 > y11)
  }
  # get all permutations of columns
  permutations <- gtools::permutations(n = ncol(mat), r = 2, 
                                       v = seq_len(ncol(mat)))
  # create a matrix to store the results
  result <- matrix(0, nrow = ncol(mat), ncol = ncol(mat))
  # apply the function to each combination
  result[permutations] <- pbapply::pbapply(permutations, 1, function(r) p2(
    mat[, r[1]], mat[, r[2]]))
  return(result)
}
