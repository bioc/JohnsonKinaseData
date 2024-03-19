\name{JohnsonKinasePWM}
\alias{JohnsonKinasePWM}
\title{Get file path to local cache of kinase PWM data}
\usage{
JohnsonKinasePWM(metadata=FALSE)
}
\arguments{
  \item{metadata}{
  When 'metadata=FALSE' the full resource is loaded, when 'metadata=TRUE' just the metadata are displayed.
  }
}
\value{
A character vector with the file path to the local cache of the PWM data.
}
\description{
Internal cache accessor function.
}
\seealso{
Use \code{\link{getKinasePWM()}} to load the data.
}
\examples{
data <- JohnsonKinasePWM()
}