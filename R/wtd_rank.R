#' Implement weighted rank correlation, as described in Jerrold H Zar, 
#' "Biostatistical Analysis" 5th ed. 
#' 
#' @param mat a matrix of data
#' @return the estimated correlation coefficients
#' @export
wtd_rank <- function(mat) {
  ranks <- apply(mat, 2, rank, ties = "average")
  # weight the ranks 
  # calculate the savage scores 
  n <- nrow(mat)
  reciprocals <- 1 / seq_len(n)
  savage <- sapply(seq_len(n), function(i) sum(reciprocals[i:n]))
  # replace each rank with the savage score 
  savages <- ranks
  savages[] <- savage[ranks]
  # calculate pearson correlation
  cor <- WGCNA::cor(savages, method = "pearson")
  return(cor)
}
