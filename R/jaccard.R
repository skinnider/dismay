#' Calculate the Jaccard index between all columns in a matrix
#' 
#' Calculate the Jaccard index between all pairs of columns in a matrix.
#' The Jaccard index (also known as the Tanimoto coefficient) is defined as the 
#' size of the intersection of two bitsets divided by the size of the union. 
#' Here, to convert a list of continuous expression values into a set of bits,
#' measurements that are missing or equal to zero are considered as zeroes and
#' all other measurements are considered as ones. 
#' 
#' @param mat a matrix of data, with samples in rows and features in columns
#' @return the Jaccard index between non-zero/missing values in each pair of 
#'   columns
#' 
#' @examples 
#' mat = matrix(c(1, rep(0, 9), rep(1, 4), rep(0, 6)), ncol = 2)
#' jaccard(mat)
#' mat = cbind(mat, c(0, rep(1, 5), rep(0, 4)))
#' 
#' @export
jaccard = function(mat) {
  present <- !is.na(mat) & mat > 0
  intersect <- crossprod(present)
  union <- nrow(mat) - crossprod(!present)
  J <- intersect / union
  return(J)
}
