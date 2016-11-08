\name{rgb2gray}
\alias{rgb2gray}
\title{ Convert RGB image to Grayscale
}
\description{
A function returns grayscale image with coefficients = c(0.3, 0.59, 0.11).
}
\usage{
rgb2gray(x, coefs=c(0.30, 0.59, 0.11))
}
\arguments{
  \item{x}{
A raster image or a matrix
}
  \item{coefs}{ R, G, B weights. Default are coefs=c(0.30, 0.59, 0.11)
}
}

\value{ A grayscale image

}
\seealso{ gray2bin
}
