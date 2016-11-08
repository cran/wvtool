\name{cryptomeria}
\alias{cryptomeria}
\docType{data}
\title{
Image Sample Dataset
}
\description{
An optical micrographs of Cryptomeria japonica
}
\usage{data("cryptomeria")}
\format{
 num [1:1079, 1:1000, 1:4] 0.886 0.89 0.863 0.639 0.424 ...
 - attr(*, "bits.per.sample")= int 8
 - attr(*, "samples.per.pixel")= int 4
 - attr(*, "sample.format")= chr "uint"
 - attr(*, "planar.config")= chr "contiguous"
 - attr(*, "compression")= chr "none"
 - attr(*, "x.resolution")= num 72
 - attr(*, "y.resolution")= num 72
 - attr(*, "resolution.unit")= chr "inch"
 - attr(*, "orientation")= chr "top.left"
 - attr(*, "artist")= chr "DP"
 - attr(*, "date.time")= chr "2016:07:29 11:38:28"
 - attr(*, "color.space")= chr "RGB"
}
\source{
Kyoto University Xylarium Database
}
\references{
http://database.rish.kyoto-u.ac.jp
}
\examples{
data(cryptomeria)
## maybe str(cryptomeria)}
