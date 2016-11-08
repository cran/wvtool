\name{glcm}
\alias{glcm}
\title{ Gray Level Co-occurence Matrix (glcm) 
}
\description{ This function supports calculating gray level co-occurence matrices from a grayscale image (< 8 bit) with requested gray level. The gray level of the source image is read from the attributes data from input TIFF file.
}
\usage{
glcm(x, t.level=4, d=1)
}
\arguments{
  \item{x}{ A gray scale image or matrix.  "x" assumes an output from readTIFF(filename, as.is=T, info=T)
}
  \item{t.level}{ A target grey level for GLCM calculation in bite. The grayscale is truncated linearly.
}
  \item{d}{ Displacement between adjacent i, j points in pixel.
}
}
\details{ The data in matrix is either inspected as images or subsequently used to calculate Haralick texture features, originally published 15 features(Haralick et al., 1973) and two additionals (Albregsten, 1995).
}
\value{ The gray level cooccurence matrices of 4 directions (theta) and their average, gray level, and displacement vector are listed.
\item{glcm }{GLCM at theta = "0", " 45",  "90",  "135" degree and "average"}
\item{level }{ numbr of gray level }
\item{d}{ length of displacement vector}
}
\references{R.M. Haralick, K. Shangmugam, Its'hak Dinstein (1973) Textural Features for Image Classification, IEEE Transactions on Systems, Man, and Cybernetics, SMC-3(6), 610-621.
	
Albregtsen F (1995) Statistical texture measures computed from gray level cooccurrence matrices. In: Technical Note, Department of Informatics, University of Oslo, Norway

K. Kobayashi, M. Akada, T. Torigoe, S. Imazu, J. Sugiyama (2015) Automated recognition of wood used in traditional Japanese sculptures by texture analysis of their low-resolution computerd tomography data. J. Wood Sci., 61, 630-640
}
\seealso{ gray2bin, rgb2gray, haralick, 
}
\examples{
data(camphora)
img <- camphora
par(mfrow=c(1,2))
lev <- 4
theta <- c(1,3)  # "th_0","th_90"
theta_c <-c("th_0","th_90")
dist <- 1
for (i in 1:2) {
	tst <- glcm(img,lev,dist)
	title <- paste(lev, "bit", " glcm ", theta_c[i], " d=", dist, sep="")
	persp(tst$glcm[[i]], theta=30, phi=30,main=title, asp=1, 
	xlab="i", ylab="j", zlab="probability")
}
}