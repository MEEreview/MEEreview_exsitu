\name{FCSex}

\alias{FCSex}

\title{FCSex}

\description{
Compute Khoury et al.'s FCSex and its components
}

\usage{
FCSex(seeds, herbarium, srange, ecoregions, bsize=50000, inrange=TRUE)
}

\arguments{
\item{seeds}{SpatVector with seed bank collection locations}
\item{herbarium}{SpatVector with herbarium collection locations}
\item{srange}{SpatRaster with the range (all values > 0 are in the range)}
\item{ecoregions}{SpatRaster with ecoregions (each ecoregion should have a unique ID)}
\item{bsize}{positive number. Buffer size in m}
\item{inrange}{logical. If \code{TRUE} buffer areas for seed locations that are outside of \code{srange} are not considered in the computation of GRSex and ERSex}
}

\note{
In \code{GapAnalysis::GRSex} the CA50 areas outside of the SDM range are considered for the seedbank area and in \code{GapAnalysis::ERSex} the number of ecoregions considered for the seedbank samples includes areas that are outside the SDM. In both cases you could end up with a score > 1. That can be corrected for, but this is inconsistent. Areas outside the range (SDM) should not be included in the analysis. They are not included in the analysis below if \code{inrange=FALSE}
}

\value{
numeric 
}

\examples{
2 + 2
}

\keyword{spatial}

