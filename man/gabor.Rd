\name{gabor.filter}
\alias{gabor.filter}
\title{
Two Dimensional Gabor Filtering in Frequency Domain
}
\description{
The function provides two dimensional Gabor function proposed by Daugman to model the spatial summation properties in the visual cortex. It returns Gabor filter in real and reciprocal space, and filtered image.}
\usage{
gabor.filter(x, lamda=5, theta=45, bw=1.5, phi=0, asp=1, disp=FALSE)
}
\arguments{
  \item{x}{A raster or matrix to be filtered }
  \item{lamda}{Wavelength of the cosine part of Gabor filter kernel in pixel. Real number greater than 2 can be used. However, lamda=2 should not be used with phase offset (phi) = -90 or 90.  }
  \item{theta}{The orientation of parallel stips of Gabor function in degree}
  \item{bw}{Half responce spatial frequency bandwidth of a Gabor filter. This relates to the ratio sigma/lamda, where sigma is the standard deviation of Gaussian factor of Gabor function. }
  \item{phi}{Phase offset of the cosine part of Gabor filter kernel in degree}
  \item{asp}{Elliptcity of the Gabor function. asp=1 means circular. For asp<1 it gives elongated parallel strip}
  \item{disp}{If this operator is TRUE, original image, gabor filter in real space domain, that in frequency domain, and filtered image will be generated }

}

\value{The function provides following four outputs.
 \item{kernel}{151x151 gabor filter kernel }
 \item{mask}{A mask with kernel in the center in a space domain}
 \item{freq_mask}{Real part of Fourier transform of the mask in spatial frequency domain}
 \item{filtered_image}{Inversed Fourier transform of FFT(img)*FFT(mask)}
}
\references{
J. G. Daugman, Uncertainty relation for resolution in space, spatial frequency, and orientation optimized by two d-dimensional visual cortical filters, J. Opt. Sppc. Am A, 2(7), 1160-1169, 1985
}
\examples{
data(cryptomeria)
img <- rgb2gray(cryptomeria)
img <- crop(img,300,300)
# simple example
test <- gabor.filter(x=img, lamda=8, theta=60, bw=1.5, phi=0, asp=0.3, disp=TRUE)
\dontrun{
# azimuthal intensity distribution with respect to the orientation of Gabor function
par(mfrow=c(2,1))
Integ <- array()
filt.img <- matrix(0, nrow(img), ncol(img))
for ( i in 1:60) {
out <- gabor.filter(x=img, lamda=8, theta=3*i, bw=1.5, phi=0, asp=0.3)
filt.img <- out$filtered_img + filt.img
Integ[i] <- sum(out$filtered_img*out$filtered_img)
}
image(rot90c(filt.img),col=gray(c(0:255)/255),asp=1,axes=FALSE, useRaster=TRUE)
x <- 1:60
plot(3*x,Integ, ty="l", ylab="integrated intensity (a.u.)", xlab="azimuthal angle (deg)")
}
}
