\name{rotate.matrix}
\alias{rotate.matrix}
\title{ Image Rotation by Bilinear Interpolation
}
\description{ Three methods to execute rotation by 1) assuming values to destination, 2) obtaining values from the source image by inverse rotation with "nearest neighbor (NN)", 3) previous procedure together with "bilinear interpolation". The default is a rotation with "bilinear Interpolation".
}
\usage{ 
rotate.matrix(x, angle, method="bilinear")
}
\arguments{
  \item{x}{
A raster image or a matrix
}
  \item{angle}{
Plus(>0) value to request clockwise rotation, while minus for anticlockwise rotation.
 }
   \item{method}{ "simple" assumes values to destination', "NN" obtains values from the source image by inverse rotation with "nearest neighbor", and "bilinear" performs  the same but with "bilinear interpolation" of the source image.
value to request clockwise rotation.
 }
}
\details{
Assuming 8-bit grayscale image as an input.
}
\value{ 
A matrix after rotation
}
\seealso{
rgb2gray
}
\examples{
data(camphora)
par(mfrow=c(2,2))
r1 <- rotate.matrix(camphora,15, method="simple")
image(rot90c(r1),asp=1,col=grey(c(0:255)/255), main= "simple", 
useRaster=TRUE, axes=FALSE)
r2 <- rotate.matrix(camphora,25, method="NN")
image(rot90c(r2),asp=1,col=grey(c(0:255)/255), main="nearest neighbour", useRaster=TRUE, axes=FALSE)
r3 <- rotate.matrix(camphora,35, method="bilinear")
image(rot90c(r3),asp=1,col=grey(c(0:255)/255),main="bilinear interpolation", 
useRaster=TRUE, axes=FALSE)
}
