\name{gatingTemplate-class}
\alias{gatingTemplate}
\alias{gatingTemplate-class}
\alias{gatingTemplate,character-method}
\alias{show}
\alias{show,boolMethod-method}
\alias{show,fcFilter-method}
\alias{show,gatingTemplate-method}
\title{a class storing the gating method and population information in a graphNEL object}
\arguments{
  \item{x}{\code{character} csv file name}

  \item{...}{other argumentss

  name: \code{character} the label of the gating template}
}
\description{
  Each cell population is stored in graph node and is
  connected with its parent population or its reference
  node for boolGate or refGate.

  It parses the csv file that specifies the gating scheme
  for a particular staining pannel.
}
\details{
  This csv must have the following columns:

  'alias': a name used label the cell population, the path
  composed by the alias and its precedent nodes (e.g.
  /root/A/B/alias) has to be uniquely identifiable.  So
  alias can not contain '/' character, which is reserved as
  path delimiter.

  'pop': population patterns of 'A+/-` or 'A+/-B+/-', which
  tells the algorithm which side (postive or negative) of
  1d gate or which quadrant of 2d gate to be kept

  'parent': the parent population alias, its path has to be
  uniquely identifiable.

  'dims': characters seperated by comma specifying the
  dimensions(1d or 2d) used for gating. It can be either
  channel name or stained marker name.

  'gating_method': the name of the gating function (e.g.
  'flowClust'). It is invoked by a wrapper function that
  has the identical function name prefixed with a dot.(e.g.
  '.flowClust')

  'gating_args': the named arguments passed to gating
  function

  'collapseDataForGating': When TRUE, data is collapsed
  (within groups if 'groupBy' specified) before gating and
  the gate is replicated across collapsed samples.  When
  set FALSE (or blank),then 'groupBy' argument is only used
  by 'preprocessing' and ignored by gating.

  'groupBy': If given, samples are split into groups by the
  unique combinations of study variable (i.e. column names
  of pData,e.g."PTID:VISITNO").  when split is numeric,
  then samples are grouped by every N samples

  'preprocessing_method': the name of the preprocessing
  function(e.g. 'prior_flowClust'). It is invoked by a
  wrapper function that has the identical function name
  prefixed with a dot.(e.g. '.prior_flowClust') the
  preprocessing results are then passed to gating wrapper
  function through 'pps_res' argument.

  'preprocessing_args': the named arguments passed to
  preprocessing function.
}
\examples{
\dontrun{
  gt <- gatingTemplate(system.file("extdata/template_tcell.csv",package = "openCyto"))
  plot(gt)
}
}

