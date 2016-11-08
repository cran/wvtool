\name{gray2bin}
\alias{gray2bin}
\title{ Conversion from Grayscale to Binary Image
}
\description{
A function provides automatic clustering-based thresholding proposed by Ohtsu, and a gray scale image is converted to binary image. Initial histogram and discriminant level of binning are displayed by his=TRUE, dis=TRUE options. A threshold value can be also set manually.
}
\usage{
gray2bin(x, auto=TRUE, th=200, his=FALSE, dis=FALSE)
}
\arguments{
  \item{x}{
A raster image or a matrix
}
  \item{auto}{
set threshold automatically (Ohtsu method) or manually
}
	\item{th}{
a threshold value used when auto=FALSE
}
  \item{his}{ 
A histogram of initial gray scale image 
}
  \item{dis}{ 
A plot of variation between classes divided by variation within classes
}
}
\value{ A requested binary image. Black is zero.
}
\references{
N. Otsu (1979) A threshold selection method from gray-level histograms, IEEE Trans. Sys., Man., Cyber., 9(1), 62-66.
}
\seealso{
rgb2gray
}
\examples{
data(camphora)
par(mfrow=c(2,3))
image(rot90c(camphora), col= gray((0:255)/255), main="camphora", asp=1, useRaster=TRUE, axes=FALSE)
out <- gray2bin(camphora, his=TRUE, dis=TRUE)
image(rot90c(out), col= gray((0:255)/255), main="binary image, auto", asp=1, 
useRaster=TRUE, axes=FALSE)
image(rot90c(gray2bin(camphora,auto=FALSE,th=100)), col= gray((0:255)/255), main="binary image, 
thresh=100", asp=1, useRaster=TRUE, axes=FALSE)
image(rot90c(gray2bin(camphora,auto=FALSE,th=180)), col= gray((0:255)/255), main="binary image, 
thresh=180", asp=1, useRaster=TRUE, axes=FALSE)
}
