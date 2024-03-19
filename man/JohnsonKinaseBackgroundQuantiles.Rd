\name{JohnsonKinaseBackgroundQuantiles}
\alias{JohnsonKinaseBackgroundQuantiles}
\title{Get file path to local cache of background score quantiles}
\usage{
JohnsonKinaseBackgroundQuantiles(metadata=FALSE)
}
\arguments{
  \item{metadata}{
  When 'metadata=FALSE' the full resource is loaded, when 'metadata=TRUE' just the metadata are displayed.
  }
}
\value{
A character vector with the file path to the local cache of background score quantiles. 
}
\description{
Internal cache accessor function.
}
\seealso{
Use \code{\link[=getScoreMaps]{getScoreMaps()}} to load the data.
}
\examples{
quants <- JohnsonKinaseBackgroundQuantiles()
}