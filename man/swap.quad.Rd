\name{swap.quad}
\alias{swap.quad}
\title{ Swapping Quadrants
}
\description{ A function maybe useful to generates power spectrum from fft output.
}
\usage{
swap.quad(x, disp=FALSE, reverse=FALSE)
}
\arguments{
  \item{x}{ output of Mod(fft(imagefile))
}
  \item{disp}{ TRUE requests to draw power spectrum 
}
  \item{reverse}{ TRUE should be used when power spectrum of N x M, where one of them is odd number. 
}
}
\value{ a matrix of power spectrum
}
\seealso{
fft, Mod
}
\examples{
data(camphora)
data(cryptomeria)
img1 <-camphora
img2 <- rgb2gray(cryptomeria)
par(mfrow=c(2,2))

image(rot90c(img1),col=gray(c(0:255)/255), main="Camphora", asp=1, 
useRaster=TRUE, axes=FALSE)
o.fft <-Mod(fft(img1))
ps <- swap.quad(o.fft)
image(rot90c(log(ps)),col=gray(c(0:255)/255), main="power spectrum", asp=1, 
useRaster=TRUE, axes=FALSE)

image(rot90c(img2),col=gray(c(0:255)/255), main="Cryptomeria", asp=1, 
useRaster=TRUE, axes=FALSE)
image(rot90c(log(swap.quad(Mod(fft(img2))))),col=gray(c(0:255)/255), 
main="power spectrum", asp=1, useRaster=TRUE, axes=FALSE)
}
