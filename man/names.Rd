\name{names}
\alias{alias}
\alias{alias,gtPopulation-method}
\alias{dims}
\alias{dims,gtMethod-method}
\alias{groupBy}
\alias{groupBy,gtMethod-method}
\alias{isCollapse}
\alias{isCollapse,gtMethod-method}
\alias{names}
\alias{names,gtMethod-method}
\alias{names,gtPopulation-method}
\alias{parameters}
\alias{parameters,gtMethod-method}
\title{get gating method name}
\arguments{
  \item{object}{\code{gtMethod}}

  \item{object}{\code{gtMethod}}

  \item{object}{\code{gtMethod}}

  \item{object}{\code{gtMethod}}

  \item{object}{\code{gtMethod}}

  \item{\code{gtPopulation}}{object}

  \item{\code{gtPopulation}}{object}
}
\value{
  \code{logical}
}
\description{
  get gating method name

  get gating method dimensions

  get parameters of the gating method/function

  When TRUE, the flow data(multiple flowFrames) is
  collapsed into one and the gating method is applied on
  the collapsed data. Once the gate is generated, it is
  thenreplicated and applied to the each single flowFrame.

  When specified, the flow data is grouped by the grouping
  variable (column names in pData). Within each group, when
  \code{isCollapse} is set to TRUE, the gating method is
  applied to the collapsed data. Otherwise, it is done
  indepentently for each indiviudal sample(flowFrame).
  Grouping variable is also used by \code{preprocessing
  method}.

  get population name

  get population alias
}
\examples{
gt <- gatingTemplate(system.file("extdata/template_tcell.csv",package = "openCyto"))

gtMthd <- getGate(gt, "/nonDebris/singlets",  "/nonDebris/singlets/lymph")
names(gtMthd)
dims(gtMthd)
parameters(gtMthd)
isCollapse(gtMthd)
groupBy(gtMthd)

gtPop <- getNodes(gt, "/nonDebris/singlets/lymph/cd3/cd4+cd8-/CD38+")
names(gtPop)
alias(gtPop)
}

