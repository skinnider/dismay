#' Fast calculation of Kendall's tau rank correlation coefficient for paired
#' non-zero observations only. Adapted from the pcaPP package source. 
#'
#' @param x A vector, a matrix or a data frame of data.
#' @param y	A vector of data.
#' @return The estimated correlation coefficient. 
cor.fk.nz <- function (x, y = NULL) {
  if (is.null(y)) {		##	x is expected to be matrix or data frame
    if (!is.matrix (x) && ! is.data.frame (x))
      stop ("x must be either numeric vector, matrix or data.frame.")
    
    p <- ncol(x)
    dn <- colnames(x)
    ret <- diag (p)
    dimnames(ret) <- list(dn, dn)

    xNA <- x
    xNA[x == 0] <- NA
    
    for (i in 1:p) {
      if (i == p)
        return(ret)
      
      ord <- order(x[, i])
      cur.x <- x[ord, i]
      for (j in (i+1):p) {
        cur.y <- x[ord, j]
        nz <- cur.x != 0 & cur.y != 0
        ret[i, j] <- ret[j, i] <- cor.fk.2d(cur.x[nz], cur.y[nz], T)
      }
    }
  } else {
    if (length (x) != length (y))
      stop ("x and y must have same length.")
    ord <- order (x)
    return (cor.fk.2d(x[ord], y[ord], cor))
  }
}
