# frSVD
Given an m by n sparse matrix A, function frsvd() can approximately computer its largest k singular values and the corresponding singular vectors. The algorithm has similar accuracy to the basic randomized SVD (rPCA) algorithm (Halko et al., 2011), but is largely optimized for sparse data. It also has good flexibility to trade off runtime against accuracy for practical usage.
