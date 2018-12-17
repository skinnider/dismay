#' Calculate the Hamming distance between all columns in a matrix
#' 
#' Calculate the Hamming distance between all pairs of columns in a matrix.
#' The Hamming distance is defined as the number of positions in two vectors
#' at which the corresponding symbols or different. Here, the Hamming distance
#' is applied to the presence/absence of each gene, whereby
#' measurements that are missing or equal to zero are considered as zeroes and
#' all other measurements are considered as ones. 
#' 
#' The matrix implementation used in \code{dismay} was obtained from:
#' \url{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
#' 
#' @param mat a matrix of data, with samples in rows and features in columns
#' @return the Hamming index between non-zero/missing values in each pair of 
#'   columns
#' 
#' @examples 
#' mat = matrix(c(1, rep(0, 9), rep(1, 4), rep(0, 6)), ncol = 2)
#' hamming(mat)
#' mat = cbind(mat, c(0, rep(1, 5), rep(0, 4)))
#' hamming(mat)
#' 
#' @export
hamming = function(mat) {
  present <- !is.na(mat) & mat > 0
  D <- t(1 - present) %*% present
  D + t(D)
}
