\name{frSVD}
\alias{frSVD}
\title{computing approximate SVD}
\description{
Given an m by n sparse matrix A, function frsvd() can approximately computer its largest k singular values and the corresponding singular vectors.
The algorithm has similar accuracy to the basic randomized SVD (rPCA) algorithm (Halko et al., 2011), but is largely optimized for sparse data. It also has good flexibility to trade off runtime against accuracy for practical usage.
}
\usage{
frSVD(A,k,q)
}
\arguments{
  \item{A}{a sparse matrix}
  \item{k}{rank k or the largest k singular value}
  \item{q}{the iteration times is equal to (q-1)/2, mostly setting 3~6}
}
\details{
nothing
}
\value{
  \item{d}{ A vector of the computed singular values.}
  \item{u}{ An m by nu matrix whose columns contain the left singular vectors.}
  \item{v}{ An n by nv matrix whose columns contain the right singular vectors.}
}
\references{
Xu Fen, Yuyang Xie, Mingye Song Wenjian Yu, Jie Tang(2018). Fast Randomized PCA for Sparse Data. Proceedings of Machine Learning Research;see also https://arxiv.org/pdf/1810.06825.pdf
}
\examples{
if(!require("Matrix")){
  install.packages("Matrix")
}
library(Matrix)
i <- sample(1:10000,2000000,replace=T)
j <- sample(1:10000,2000000,replace=T)
x <- rep(1,2000000)
A <- sparseMatrix(i, j, x = x)
frSVD(A,100,3)
}
