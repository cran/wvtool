\name{dec2bin}
\alias{dec2bin}
\title{
	Decimal to Binary Convertion 
}
\description{ 
	A function returns binary number from decimal number. The function is internally used in lbp function.
}
\usage{
	dec2bin(x, digit=8)
}
\arguments{
  \item{x}{
	A decimal integer.
}
  \item{digit}{
	A length of binary sequence.
}
}
\value{ 
	A binary number in array.
}
\seealso{ bin2dec, lbp
}