\name{ppMethod-class}
\alias{ppMethod-class}
\title{A class to represent a preprocessing method.}
\description{
  It extends \code{gtMethod} class.
}
\examples{
\dontrun{
     gt <- gatingTemplate(system.file("extdata/tcell.csv",package = "openCyto"))
     ppMethod(gt, '3', '4')
}
}

