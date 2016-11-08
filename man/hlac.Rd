\name{hlac}
\alias{hlac}
\title{ Higher Order Local Autocorrelation (HLAC)
}
\description{Feature extraction for practical vision system, whose features are shift-invariant and additive. The function gives zero to the eighth order cases, represented by 223 mask patterns of 3 x 3 within a 2r+1 x 2r+1 (r >=1) displacement region. 
}
\usage{
hlac(x, r=1, disp=FALSE)
}
\arguments{
  \item{x}{A binary or gray image or matrix
}
  \item{r}{Displacement vector r for 3 x 3 mask pattern
}
  \item{disp}{If TRUE, function saves 223 filtered images in one matrix.
}
}
\details{ The feature parameter should be a list. The function returns 1, 4, 20, 45, 62, 54, 28, 8,1 features and corresponding filtered images if disp is TRUE.
}
\value{HLAC features or the corresponding image with requested HLAC measures.
\item{features }{Numerical output of 0 to 8th order masks}
\item{mat }{A large matrix of 223 images expanded in a row }
}
\references{ N.Otsu, T. Kurita (1988) A New Scheme for Practical Flexible and Intelligent Vision Systems, In: Proc. Machine Vision Application(MVA), 431-435.
}
\seealso{rgb2gray, gray2bin, glcm, lbp
}
\examples{
# features plot and the corresponding image presentation
data(camphora)
tmp <- hlac(gray2bin(camphora), 2, disp=TRUE)
par(mfrow=c(2,2))
plot(unlist(tmp$features), main="HLAC histogram")
image(rot90c(matrix(tmp$mat[2,],tmp$row,tmp$col)),
	col  = gray((255:0)/255), main="2", useRaster=TRUE, asp=1, axes=FALSE)
image(rot90c(matrix(tmp$mat[23,],tmp$row,tmp$col)),
	col  = gray((255:0)/255), main="23", useRaster=TRUE, asp=1, axes=FALSE)
image(rot90c(matrix(tmp$mat[156,],tmp$row,tmp$col)),
	col  = gray((255:0)/255), main="156", useRaster=TRUE, asp=1, axes=FALSE)
}