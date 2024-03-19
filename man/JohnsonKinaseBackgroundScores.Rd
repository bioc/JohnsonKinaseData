\name{JohnsonKinaseBackgroundScores}
\alias{JohnsonKinaseBackgroundScores}
\title{Get file path to local cache of background scores}
\usage{
JohnsonKinaseBackgroundScores(metadata=FALSE)
}
\arguments{
  \item{metadata}{
  When 'metadata=FALSE' the full resource is loaded, when 'metadata=TRUE' just the metadata are displayed.
  }
}
\value{
A character vector with the file path to the local cache of the background scores. 
}
\description{
Internal cache accessor function.
}
\seealso{
Use \code{\link{getBackgroundScores()}} to load the data.
}
\examples{
scores <- JohnsonKinaseBackgroundScores()
}