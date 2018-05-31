CytoGate
========

Plugin Gating Functions for openCyto

Installation
============

**CytoGate** can be installed from Github

GitHub
------

``` r
library(devtools)
install_github("DillonHammill/CytoGate")
```

Plugin Example
==============

To demonstrate the use of **CytoGate** the following example uses the **gate\_draw** plugin for openCyto to draw polygonGates manually using the mouse.

Register openCyto Plugin
------------------------

``` r
library(openCyto)
library(CytoGate)

registerPlugins(fun = gate_draw, methodName = "DrawGate")
listgtMethods()
```

Construct GatingSet
-------------------

``` r
data(fs, package = "CytoGate")
gs <- GatingSet(fs)
```

Construct GatingTemplate Entry
------------------------------

``` r
library(flowDensity)

template <- add_pop(
gs, alias = "Lymphocytes", pop = "+", parent = "root", dims = "FSC-A,SSC-A", 
gating_method = "DrawGate", gating_args = "subSample=20000", collapseDataForGating = TRUE,
groupBy = 2
)
```

    ## Draw 2D polygon gate around population.

Plot Gating Result
------------------

``` r
library(ggcyto)
ggcyto(gs, subset = "root", aes(x = "FSC-A", y = "SSC-A")) + geom_hex(bins = 100) + geom_gate("Lymphocytes") + geom_stats()
![](https://github.com/DillonHammill/CytoGate/tree/master/README_files/figure-markdown_github-ascii_identifiers)
```

Dillon Hammill, BMedSci (Hons) <br /> Ph.D. Scholar <br /> The Parish Group – Cancer & Vascular Biology <br /> ACRF Department of Cancer Biology and Therapeutics <br /> The John Curtin School of Medical Research <br /> ANU College of Medicine, Biology and the Environment <br /> The Australian National University <br /> Acton ACT 2601 <br /> <Dillon.Hammill@anu.edu.au>
