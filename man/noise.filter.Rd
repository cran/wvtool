\name{noise.filter}
\alias{noise.filter}
\title{
Median, Mean and Gaussian Filter
}
\description{
A funtion provides three kinds of noise reduction on an image, "median", "mean", and "gaussian". A typical pre-processing step to improve the results of later processing for example, glcm-haralick analysis.
}
\usage{
noise.filter(x, n=3, method="median")
}
\arguments{
  \item{x}{
A raster image or a matrix
}
  \item{n}{
filter size is given by n x n. Default is 3 x 3. Number has to be an odd number. For gaussian filter, only 3 or 5 is available.
}
  \item{method}{
"median", "mean", and "gaussian" can be selected. Default is "median".
}
}
\value{ A raster or a matrix
}
\references{
T.S. Huang, G.J. Yang, G.Y. Tang (1979) A fast two-dimensional median filtering algorithm, IEEE transactions, Acoustics, Speech and Signal Processing, 27, 13-18.

}
\seealso{
glcm
}
\examples{
data(camphora)
camphora <- crop(camphora,200,200)
par(mfrow=c(2,2))
image(rot90c(noise.filter(camphora,3,"median")),col=gray(c(0:255)/255), 
main="median", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(noise.filter(camphora,3,"mean")),col=gray(c(0:255)/255), 
main="mean", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(noise.filter(camphora,3,"gaussian")),col=gray(c(0:255)/255), 
main="gaussian", useRaster=TRUE, axes=FALSE, asp=1)
}
