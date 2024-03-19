\name{JohnsonKinaseAnnotation}
\alias{JohnsonKinaseAnnotation}
\title{Get file path to local cache of annotation data for kinase PWMs}
\usage{
JohnsonKinaseAnnotation(metadata=FALSE)
}
\arguments{
  \item{metadata}{
  When 'metadata=FALSE' the full resource is loaded, when 'metadata=TRUE' just the metadata are displayed.
  }
}
\value{
A character vector with the file path to the local cache of the annotation data.
}
\description{
Internal accessor function.
}
\seealso{
Use \code{\link{getKinaseAnnotation()}} to load the data.
}
\examples{
anno <- JohnsonKinaseAnnotation()
}