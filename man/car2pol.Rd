\name{car2pol}
\alias{car2pol}
\title{
	Polar Transformer -Cartesian to Polar Coordinates-
}
\description{
The function converts images to polar coordintates. The polar transformation is useful for unwarpping images which have a generally round object. From power specturm for example, one may generate radial integration profile or azimuthal intensity distribution. Default is "bilinear" interpolation.
}
\usage{
car2pol(x, method="bilinear")
}
\arguments{
\item{x}{ 
A raster image or a matrix.
}
\item{method}{ "NN" Nearest neighbour method, "bilinear" Bilinear interpolation.
}
}
\value{ A matrix in polar coordinate system of the requested image 
\item{pol.img}{Radial distance corresponds to the shorter side of requested image, and polar angle covers 0 to 360 degrees.}
}
\seealso{
integ.profile
}
\examples{
data(camphora)
par(mfrow=c(1,2))
image(rot90c(camphora), col=gray(c(0:255)/255), main="camphora, original", 
	    useRaster=TRUE, asp=1, axes=FALSE)
img <- car2pol(camphora, method="bilinear")
image(rot90c(img), col=gray(c(0:255)/255), main="camphora, polar (bilinear)", 
xlab="radial distance(pixel)",ylab="angle(deg)", useRaster=TRUE, asp=1, axes=FALSE)
}
