\name{integ.profile}
\alias{integ.profile}
\title{
Simple Integration for Making Profile
}
\description{
A function returns integrated line profile. It crops rectangular area from a requested size and project and integrate pixel values either to horizontal or vertical axis. 
When used with a matrix in polar coordinate (car2pol) calculated from power spectrum (power.spec) of an image, the function provides radial integration or azimuthal integration that are useful for diffraction analysis.
}
\usage{
integ.profile(x, axis="H", h=c(20, 50), v=c(30, 120), disp=FALSE)
}
\arguments{
  \item{x}{
A raster image or a matrix
}
  \item{axis}{ Axis to project. 
H : Projection to horizontal axis or radial distance  ( in polar coordinate).
V : Projection to vertical axis or azimuthal angle ( in polar coordinate).
}
  \item{h}{ c(h1,h2) :  A horizontal or radial (in polar coordinate) range for integration.
  	}
  \item{v}{ c(v1,v2) :  A vertical or azimuthal (in polar coordinate) range for integration.
}
  \item{disp}{
Plot calculated profile. Default is FALSE.
}
}
\details{
The row and column corresponds to horizontal and vertical axes, respectively.
}
\value{ An array of requested line profile
}
\seealso{
swap.quad, car2pol, fft, Mod
}
\examples{
data("camphora")
img <- camphora
par(mfrow=c(2,2))
image(rot90c(img),col=gray(c(0:255)/255), useRaster=TRUE, main="camphora",asp=1, axes=FALSE)
integ.profile(img, axis="H", h=c(1,nrow(img)) , v=c(1,ncol(img)), disp=TRUE)
integ.profile(img, axis="V", h=c(1,nrow(img)) , v=c(1,ncol(img)), disp=TRUE)
ps <- log(swap.quad(Mod(fft(img))))
pol <- car2pol(ps)
image(rot90c(ps),col=gray(c(0:255)/255), useRaster=TRUE,main="power spectrum",asp=1, axes=FALSE)
image(rot90c(pol),col=gray(c(0:255)/255), useRaster=TRUE,main="polar map",asp=1)
integ.profile(pol, axis="H", h=c(10,200) , v=c(0,90), disp=TRUE)
integ.profile(pol, axis="V", h=c(70,100) , v=c(0,360), disp=TRUE)
}
