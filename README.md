openCyto [<img src="logo_mid.png"/>](http://github.com/RGLab/openCyto)
========


An R package that providing an automated data analysis pipeline for flow cytometry.

## Documentation Resources
The [github-pages](http://opencyto.org) site documentation is a bit outdated.

The [package vignettes](http://www.bioconductor.org/packages/devel/bioc/html/openCyto.html) are the current best resource to 
get started with the openCyto.

- [An introduction to the openCyto package](https://bioconductor.org/packages/devel/bioc/vignettes/openCyto/inst/doc/openCytoVignette.html)
- [How to use different auto gating functions](https://bioconductor.org/packages/devel/bioc/vignettes/openCyto/inst/doc/HowToAutoGating.html)
- [How to write a csv gating template](https://bioconductor.org/packages/devel/bioc/vignettes/openCyto/inst/doc/HowToWriteCSVTemplate.html)

The `add_pop()` API is a good interactive approach to building up a template, population by population. It takes arguments found in the `csv` template, performs the gating of a `GatingSet` or `GatingHierarchy` and returns a line of text that can be added to a `csv` template.

