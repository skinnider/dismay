cor.fk.2d <- function (x, y, cor) {
  if (length (x) != length (y))
    stop ("x and y must have same length.")
  
  ret <- .C(pcaPP:::C_kendallNlogN, PACKAGE = "pcaPP", NAOK = FALSE, DUP = TRUE	
            ##	20130322 set DUP = TRUE - kendallNlogN implementation modifies x & y vectors!!
            , as.double (x)
            , as.double (y)
            , as.integer (c (length (x), cor))
            , ret = double (1)
  )
  return (ret$ret)
}
