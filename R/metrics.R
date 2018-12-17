#' Get all measures of distance or similarity implemented in \code{dismay}
#' 
#' @return all measures of distance or similarity implemented in \code{dismay}
#' @export
metrics = function() {
  c('pearson', 'spearman', 'kendall', 'bicor', 'binomial', 
    'cosine', 'jaccard', 'canberra', 'euclidean', 'manhattan', 'weighted_rank',
    'hamming', 'dice', 'phi_s', 'rho_p', 'zi_kendall', 'MI')
}

#' Get metrics bounded by the range [-1, 1]
#' 
#' Test whether the values for a metric of interest lie within the range
#' [-1, 1]. This returns true for Pearson/Spearman/Kendall correlations, 
#' biweight midcorrelation, zero-inflated Kendall's tau, mutual information
#' (as implemented in WGCNA), the \code{rho_p} measure of proportionality, 
#' weighted rank correlation, and cosine similarity.  
#' 
#' @return true if the metric is within the range [-1, 1], false otherwise
#' @export
is_bounded = function(metric) {
  metric %in% c("pearson", "spearman", "kendall", "bicor", "zi_kendall",
                "MI", "rho_p", "weighted_rank", "cosine")
}