\name{ex_cvs}

\alias{ex_cvs}

  
\title{ex-situ conservation score}

\description{
get the ex-situ conservation score for a taxon
}

\usage{
ex_cvs(regions, sample, env=NULL, adjust=TRUE, minssize=10)
}

\arguments{
\item{regions}{SpatVector of polygons}
\item{sample}{SpatVector of points}
\item{env}{SpatRaster with environmental data or NULL}
\item{adjust}{logical. If \code{TRUE} then the cover is used?}
\item{minssize}{numeric}
}

\value{
numeric 
}

\examples{
2 + 2
}

\keyword{spatial}

