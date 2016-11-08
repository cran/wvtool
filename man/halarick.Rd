\name{haralick}
\alias{haralick}
\title{ Haralick Texture Features Calculated from GLCM
}
\description{ A function returns 15 Halarick features for 4 directions, their average and range. 
}
\usage{
haralick(x)
}
\arguments{
  \item{x}{
output of glcm() function from a TIFF data}
}
\details{ 15 outputs are
#1 Angular Second Moment / Homogeniety "asm"	
#2 Contrast "con"
#3 inverse Difference Moment "idm"
#4 Entropy "ent"
#5 Correlation 	"cor"
#6 Variance in Haralick 1973 "var"
#7 Sum Average "sav"
#8 Sum Entropy "sen"
#9 Difference Entropy "den"
#10 Difference Variance "dva"
#11 Sum Variance "sva"
#12 Information Measures of Correlation "f12" 
#13 Information Measures of Correlation "f13"
#14 Cluster Shade "sha"
#15 Cluster prominence "pro", respectively

}
\value{ A matrix of angles and features
}
\references{R.M. Haralick, K. Shangmugam, Its'hak Dinstein (1973) Textural Features for Image Classification, IEEE Transactions on Systems, Man, and Cybernetics, SMC-3(6), 610-621.
	
Albregtsen F (1995) Statistical texture measures computed from gray level cooccurrence matrices. In: Technical Note, Department of Informatics, University of Oslo, Norway

K. Kobayashi, M. Akada, T. Torigoe, S. Imazu, J. Sugiyama (2015) Automated recognition of wood used in traditional Japanese sculptures by texture analysis of their low-resolution computerd comography data. J. Wood Sci., 61, 630-640.
}
\seealso{
glcm
}
\examples{
data(camphora)
haralick(glcm(camphora,6,1))
}
