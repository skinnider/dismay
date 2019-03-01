# README

The `dismay` function provides a common interface to calculate distance metrics or measures of association for matrices. Implemented metrics include standard correlations (Pearson, Spearman, Kendall), weighted rank correlation, biweight midcorrelation, a zero-inflated estimator of Kendall's tau, mutual information, the gene co-dependency index, Jaccard distance, cosine similarity, the Sorenson-Dice coefficient, Hamming distance, Canberra, Euclidean, and Manhattan distances, and two symmetric measures of proportionality. 

## System requirements

This code requires R (>= 3.4.0) and the following dependencies:

``` 
	arules (>=1.6-1),
	gtools (>= 3.5.0),
	lsa (>= 0.73.1),
	Matrix (>= 1.2.10),
	pbapply (>= 1.3.4),
	pcaPP (>= 1.9.72),
	philentropy (>= 0.0.3),
	propr (>= 3.1.8),
	reshape2 (>= 1.4.2),
	WGCNA (>= 1.51),
	Rdpack
```

The main function takes as input a numeric matrix with M columns and returns a similarity matrix of dimensions MxM. Further details are provided in the documentation of `dismay.R`. 

## Installation

To install `dismay` on your system, download the directory and install with devtools:

```
devtools::install("/path/to/dismay")
```

Or, install directly from GitHub:

```
devtools::install_github("skinnider/dismay")
```

This should take no more than a few minutes.

## Demonstration

Generating simulated gene coexpression networks (or cell-cell similarity matrices) with different measures of association using random data is easy, with a single interface to all seventeen measures of association:

```
mat = matrix(rnorm(100), ncol = 10, dimnames = list(paste0("cell_", 1:10), 
   paste0("gene_", letters[1:10])))
mat[mat < 0] = 0
tc = dismay(mat, 'jaccard')
cos = dismay(mat, 'cosine')
```
