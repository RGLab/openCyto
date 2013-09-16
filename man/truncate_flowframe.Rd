\name{truncate_flowframe}
\alias{truncate_flowframe}
\title{Removes any observation from the given flowFrame object that has values
outside the given range for the specified channels}
\usage{
  truncate_flowframe(flow_frame, channels, min = NULL,
    max = NULL)
}
\arguments{
  \item{flow_frame}{a \code{flowFrame} object}

  \item{min}{a numeric vector that sets the lower bounds
  for data filtering}

  \item{max}{a numeric vector that sets the upper bounds
  for data filtering}

  \item{channels}{\code{character} specifying which channel
  to operate on}
}
\value{
  a \code{flowFrame} object
}
\description{
  The minimum/maximum values are ignored if \code{NULL}.
}
\examples{
library(flowClust)
data(rituximab)
# Consider the range of values for FSC.H and SSC.H
summary(rituximab)

# Truncates any observations with FSC.H outside [100, 950]
rituximab2 <- truncate_flowframe(rituximab, channels = "FSC.H", min = 100,
max = 950)
summary(rituximab2)
# Next, truncates any observations with SSC.H outside [50, 1000]
rituximab3 <- truncate_flowframe(rituximab2, channels = "SSC.H", min = 50,
max = 1000)
summary(rituximab3)

# Instead, truncates both channels at the same time
rituximab4 <- truncate_flowframe(rituximab, channels = c("FSC.H", "SSC.H"),
min = c(100, 50), max = c(950, 1000))
summary(rituximab4)
}

