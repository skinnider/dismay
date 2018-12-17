#' Calculate the Dice distance between all columns in a matrix
#' 
#' Calculate the Dice distance between all pairs of columns in a matrix.
#' The Dice distance is defined as the number of positions in two vectors
#' at which the corresponding symbols or different. Here, the Dice distance
#' is applied to the presence/absence of each gene, whereby
#' measurements that are missing or equal to zero are considered as zeroes and
#' all other measurements are considered as ones. This is a wrapper for the
#' philentropy distance function, that converts a numeric matrix to a 
#' presence/absence matrix.
#' 
#' @param mat a matrix of data, with samples in rows and features in columns
#' @return the Dice coefficient between non-zero/missing values in each pair of 
#'   columns
#' 
#' @examples 
#' mat = matrix(c(1, rep(0, 9), rep(1, 4), rep(0, 6)), ncol = 2)
#' dice(mat)
#' mat = cbind(mat, c(0, rep(1, 5), rep(0, 4)))
#' dice(mat)
#' 
#' @export
dice = function(mat) {
  present = !is.na(mat) & mat > 0
  as.matrix(arules::dissimilarity(t(present), method = 'dice'))
}
