\name{registerGatingFunction}
\alias{registerGatingFunction}
\title{Register a gating function with OpenCyto}
\usage{
  registerGatingFunction(fun = NA,
    methodName = "myGatingMethod", dep = NA)
}
\arguments{
  \item{fun}{\code{function} to be registered}

  \item{methodName}{\code{character} name of the gating
  method}

  \item{dep}{\code{character} name of the library
  dependency required for the plugin gating method to
  work.}
}
\value{
  \code{logical} TRUE if successful and prints a message.
  FALSE otherwise.
}
\description{
  Function registers a new gating method with openCyto so
  that it may be used in the csv template.
}
\details{
  The \code{fun} argument should be a wrapper function
  definition for the gating method, with formal arguments
  fr a flowFrame, pp_res a pre-processing result, xChannel
  (optional character vector), yChannel (required character
  vector), and a filter_id (character), and an ellipsis ...
  for additional parameters. The function must return a
  filter (i.e. polygonGate or other instance) from
  flowCore. Not all formal parameters need to be used.
  Additional arguments are passed via the ... and can be
  processed in the wrapper
}

