\name{cc.label}
\alias{cc.label}
\title{
Connected Component Labelling.
}
\description{
A function labels the connected components in a binary image. For example, it can be used for statistical analysis of tracheids (see examples).
}
\usage{
cc.label(x, connect=8, inv=FALSE, img.show=FALSE,text.size=0.3)
}
\arguments{
  \item{x}{
A binary image (A matrix with 0 and 1)
}
  \item{connect}{
8-connectivity or 4-connectivity. Default is 8-connectivity.
}
  \item{inv}{
inverse the binary image x before labelling. Labelling the connected area with 0 when this is TRUE.
}
  \item{img.show}{
If this is TRUE, the image with labelling numbers are shown.
}
  \item{text.size}{
the size of labelling numbers used when img.show=TRUE.
}
}
\details{
Labelling the connected components with pixels equal to 1 (white) in a binary image (If pixels equal to 0 (black) should be labelled, select inv=TRUE). The function returns the labelled image and the statistical data of the labelled components.
}
\value{
a list with 2 components (a matrix and a dataframe)
	\item{image}{
A matrix with labels
}
	\item{summary}{
A dataframe summarizing the labelled components with 8 following variables. 
}
	\item{summary$label}{labelling numbers, area: area of each component
}
	\item{summary$aveX, summary$aveY}{center position of each component
}
	\item{summary$dX, summary$dY}{width and height of each component
}
	\item{summary$edge}{If the component is on the edge of the image, this value is 1, otherwise 0.
}
}
\examples{
\dontrun{
data(cryptomeria)
img <- rgb2gray(cryptomeria)
img.c <- crop(img,300,300)
img.bin <- gray2bin(img.c, auto=FALSE, th=180)
par(mfrow=c(2,2))
test <- cc.label(img.bin, connect=8, img.show=TRUE)
hist(test$summary$area,main="histogram of area")
hist(test$summary$dX,main="histogram of dX")
hist(test$summary$dY,main="histogram of dY")
}
}