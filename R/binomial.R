#' Calculate binomial distribution P-values between all columns in a matrix
#' 
#' Calculate P-values from the binomial distribution and take their negative
#' logarithms as an indicator of coexpression, as proposed by Mohammadi et al.
#' 
#' @param mat a matrix of data, with samples in rows and features in columns
#' @return the binomial P-values between non-zero/missing values in each pair of 
#'   columns
#' 
#' @export
binomial = function(mat) {
  present = !is.na(mat) & mat > 0
  size = nrow(mat)
  q = crossprod(present)
  prob = crossprod(t(colMeans(present)))
  cor = pbinom(q, size, prob, lower.tail = F)
  dimnames(cor) = list(colnames(mat), colnames(mat))
  -log10(cor)
}
