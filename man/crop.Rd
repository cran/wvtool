\name{crop}
\alias{crop}
\title{
Image cropping
}
\description{
image cropping from the center.
}
\usage{
crop(x, width=300, height=300, shift=c(0,0))
}
\arguments{
  \item{x}{
A raster or a matrix}
  \item{width}{
width for cropping}
  \item{height}{
height for cropping}
  \item{shift}{
shift of the cropped position from the center}
 }
\value{ A raster or a matrix
}
\examples{
data(camphora)	
par(mfrow=c(2,2))
image(rot90c(camphora),col=gray(c(0:255)/255), main="original", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(crop(camphora,200,100)),col=gray(c(0:255)/255),
main="cropped from the center", useRaster=TRUE, axes=FALSE, asp=0.5)
image(rot90c(crop(camphora,200,200)),col=gray(c(0:255)/255),
main="cropped from the center", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(crop(camphora,200,200,shift=c(50,50))),col=gray(c(0:255)/255),
main="cropped from shifted position c(50,50)", useRaster=TRUE, axes=FALSE, asp=1)
}
