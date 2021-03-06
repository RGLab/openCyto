% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{gtMethod-class}
\alias{gtMethod-class}
\alias{gtMethod}
\title{A class to represent a gating method.}
\description{
A gating method object contains the specifics for generating the gates.
}
\section{Slots}{

 \describe{

     \item{name}{ a \code{character} specifying the name of the gating method}
     \item{dims}{ a \code{character} vector specifying the dimensions (channels or markers) of the gate}

     \item{args}{ a \code{list} specifying the arguments passed to gating function}

     \item{groupBy}{ a \code{character} or \code{integer} specifying how to group the data.
                  If \code{character}, group the data by the study variables (columns in \code{pData}).
                   If \code{integer},  group the data by every \code{N} samples.
                  }

     \item{collapse}{ a \code{logical} specifying wether to collapse the data within group before gating.
                 it is only valid when \code{groupBy} is specified}
 }
}

\examples{
 \dontrun{
     gt <- gatingTemplate(system.file("extdata/gating_template/tcell.csv",package = "openCyto"))
     gh_pop_get_gate(gt, '2', '3')
}
}
