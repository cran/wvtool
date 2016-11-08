\name{camphora}
\alias{camphora}
\docType{data}
\title{
Image Sample Dataset
}
\description{
An optical micrographs of Cinnamomum camphora
}
\usage{data("camphora")}
\format{
  The format is:
 num [1:486, 1:518] 0.275 0.337 0.765 0.937 0.933 ...
 - attr(*, "bits.per.sample")= int 8
 - attr(*, "samples.per.pixel")= int 1
 - attr(*, "sample.format")= chr "uint"
 - attr(*, "planar.config")= chr "contiguous"
 - attr(*, "compression")= chr "none"
 - attr(*, "x.resolution")= num 32
 - attr(*, "y.resolution")= num 32
 - attr(*, "resolution.unit")= chr "inch"
 - attr(*, "orientation")= chr "top.left"
 - attr(*, "color.space")= chr "black is zero"}
\source{
Kyoto University Xylarium Database
}
\references{
http://database.rish.kyoto-u.ac.jp
}
\examples{
data(camphora)
## maybe str(camphora) 
}
