\name{edge.detect}
\alias{edge.detect}
\title{
Canny and Sobel Edge detector.
}
\description{
A function detects edges in images by Canny or Sobel operator. Sobel provides approximate intensity of gradients for each pixels, while Canny provides a binary image with thin edges.
}
\usage{
edge.detect(x, thresh1=1, thresh2=15, noise="gaussian", noise.s=3, method="Canny")
}
\arguments{
  \item{x}{
A raster image or a matrix
}
  \item{thresh1}{
low threshold for edge tracking by hysteresis (0-100). Only used for "Canny" edge detector.
}
  \item{thresh2}{
high threshold for edge tracking by hysteresis (0-100). Only used for "Canny" edge detector.
}
  \item{noise}{
a method for noise reduction. "gaussian", "median", and "mean" filters are available. Default is "gaussian".
}
  \item{noise.s}{
filter size for noise reduction (3 or 5). Default is 3.
}
  \item{method}{
"Canny" and "Sobel" can be selected. Default is "Canny".
}

}
\details{
Canny edge detector has four steps.
1. noise reduction/
2. finding the gradient in images by Sobel operator/
3. Non-maximum suppression/
4. Hysteresis threshold.
When the method "Sobel" is selected, only step 1 and 2 will be done.
}
\value{
A raster or a matrix}
\seealso{
noise.filter
}
\examples{
\dontrun{
data(camphora)	
data(cryptomeria)
cryptomeria <- rgb2gray(cryptomeria)
img.c1 <- crop(camphora,200,200)
img.c2 <- crop(cryptomeria,300,300)
par(mfrow=c(2,2))
image(rot90c(edge.detect(img.c1,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
 method="Canny")),col=gray(c(0:255)/255), main="Canny", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(edge.detect(img.c1,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
 method="Sobel")),col=gray(c(0:255)/255), main="Sobel", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(edge.detect(img.c2,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
 method="Canny")),col=gray(c(0:255)/255), main="Canny", useRaster=TRUE, axes=FALSE, asp=1)
image(rot90c(edge.detect(img.c2,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
 method="Sobel")),col=gray(c(0:255)/255), main="Sobel", useRaster=TRUE, axes=FALSE, asp=1)
}
}