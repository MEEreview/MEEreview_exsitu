\name{adjust_range}

\alias{adjust_range}

  
\title{adjust_range}

\description{
Adjust the range
}

\usage{
adjust_range(x, sp, land, CAmin=50000, CAmax=250000)
}

\arguments{
\item{x}{SpatRaster with species range predicted with a SDM (NA where absent)}
\item{sp}{SpatVector with species occurrences}
\item{land}{SpatVector or SpatVector with land areas (or sea areas for marine species) to mask out areas in the sea (or on land)}
\item{CAmin}{positive number. Buffer size around each occurrence, to be included into the range}
\item{CAmax}{positive number. Buffer size around each occurrence, to be excluded from the range}
}

\value{
numeric 
}

\examples{
2 + 2
#adjust_range( ) 
}

\keyword{spatial}

