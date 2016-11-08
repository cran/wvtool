\name{lbp}
\alias{lbp}
\title{ Local Binary Patterns (lbp)
}
\description{ Calculate local binary patterns from a grayscale image
}
\usage{
lbp(x, r=1)
}
\arguments{
  \item{x}{ A raster image or a matrix
}
  \item{r}{ displacement vector in 8 direction. r=1 means c(-1, 0, -1,1, 0,1,1,1,1,0,1,-1,0,-1,-1,-1) r=2 means  c(-2, 1, -1,2, 1,2,2,1,2,-1,1,-2,-1,-2,-2,-1)
}
}
\details{ The LBP operator was originally designed for texture description. The operator assigns a label to every pixel of an image by thresholding the 3x3-neighborhood of each pixel with the center pixel value and considering the result as a binary number (gives 0 if each pixel is smaller than the center, otherwise 1). Then, the histogram of the labels can be used as a texture descriptor. The circular (8,r=1), and (8,r=2) neighborhoods are considered. The function assumes 8-bit grayscale image as an input.
}
\value{
\item{lbp.u2}{a matrix (image) returned from requested LBP u2 operation}
\item{lbp.ori}{a matrix (image) returned from requested LBP operation}
}
\references{
	T. Ojala, M. Pietikainen, and D. Harwood (1994), Performance evaluation of texture measures with classification based on Kullback discrimination of distributions, Proceedings of the 12th IAPR International Conference on Pattern Recognition(ICPR 1994), vol.1, pp. 582-585.
	
	T. Ojala, M. Pietikainen, and D. Harwood (1996) A Comparative Study of Texture Measures with Classification Based on Feature Distributions, Pattern Recognition, vol. 29, no. 1, 51-59.
	
}
\note{ 
A local binary pattern is called uniform if the binary pattern contains at most two 0-1 or 1-0 transitions. In calculation, the histogram has 58 separate bins for uniform patterns, and all other non-uniform patterns are assigned to one single bin. Thus, the length of the features reduces from 256 to 59. This function returned both the lbp patterns before (lbp.ori) and after (lbp.u2) considering the uniform patterns.
}
\seealso{ rgb2bin, hlac
}
\examples{
\dontrun{
data(camphora)
par(mfrow=c(2,2))
r1 <- lbp(camphora,1)
image(rot90c(r1$lbp.u2),col = gray((0:58)/58), main="lbp.u2 (r=1, 8 points)", useRaster=TRUE, 
asp=1, axes=FALSE)
image(rot90c(r1$lbp.ori),col = gray((0:255)/255), main="lbp.ori (r=1, 8 points)", useRaster=TRUE, 
asp=1, axes=FALSE)
hist(r1$lbp.u2,breaks=59, main="Histogram of lbp.u2")
hist(r1$lbp.ori,breaks=256, main="Histogram of lbp.ori")
}
}