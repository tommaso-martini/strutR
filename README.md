# strutR — Struthers Wetting-Front Model for Infiltration, Redistribution, and Drainage

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxx)

<!-- badges: start -->
[![R-CMD-check](https://github.com/tommaso-martini/strutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tommaso-martini/strutR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

[![pkgdown](https://github.com/tommaso-martini/strutR/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/tommaso-martini/strutR/actions/workflows/pkgdown.yaml)

The `strutR` package provides a complete, transparent, and open-source R implementation of the multiple wetting-front gravitational infiltration and redistribution model originally introduced by:

Struthers, I., Hinz, C., & Sivapalan, M. (2006). *A multiple wetting front gravitational infiltration and redistribution model for water balance applications.* Water Resources Research, 42(6), W06406. https://doi.org/10.1029/2005WR004645

---

## Overview

The Struthers model describes infiltration, redistribution, and drainage in a one-dimensional unsaturated soil column through the dynamics of moving wetting fronts. Each front represents a discontinuity in soil water content, moving downward as infiltration proceeds and merging when moisture contrasts become small. The model assumes gravity-driven flow between fronts, with conductivity governed by a Brooks–Corey law and water mass conserved throughout redistribution and drainage.

This R implementation, `strutR`, is designed for reproducibility, transparency, and compatibility with stochastic and deterministic hydrological modelling. It can be used both as a standalone infiltration–redistribution solver and as a lower boundary module coupled to surface or root-zone models (for example, a two-bucket system).

Main features include:

- Accurate representation of infiltration–redistribution–drainage (IRD) processes based on the Struthers et al. (2006) formulation.
- Explicit handling of multiple fronts with dynamic merging and drainage to the water table or lower boundary.
- Brooks–Corey conductivity formulation:

  $$ K(\theta) = K_s \left(\frac{\theta - \theta_r}{\theta_s - \theta_r}\right)^{1/\beta}. $$
  
- Modular structure allowing easy coupling with evapotranspiration and root-zone models.
- Full control over gravitational and capillary terms, numerical tolerances, and merging behaviour.
- Vignettes and utilities for plotting infiltration profiles, tracking front depths, and computing depth-averaged soil moisture.

This implementation aims to provide a faithful and research-ready translation of the original conceptual model, suitable for hydrological process analysis, ecohydrology, soil physics education, and stochastic soil moisture studies.

---

## Installation

### From GitHub

```r
# install.packages("devtools")  # if not already installed
devtools::install_github("tommaso-martini/strutR")
library(strutR)
```

### From a local source folder

```r
devtools::install_local(".")
library(strutR)
```

After installation, all core functions and vignettes can be accessed using `help(package = "strutR")`.

---

## Quick Example

A minimal infiltration–redistribution–drainage (IRD) simulation using the Struthers solver:

```r
# Initial conditions
fronts <- list(list(theta = 0.10, x = 2.0))

# Parameters
theta_r <- 0.05
theta_s <- 0.45
beta    <- 1/8
Ks      <- 1.0
L       <- 2.0
delta   <- 0.0
dt      <- 1/24
f_t     <- 0.01  # 10 mm/day

# Run a single redistribution step
res <- struthers_redistr_under(
  fronts = fronts,
  theta_r = theta_r, theta_s = theta_s,
  beta = beta, Ks = Ks,
  L = L, delta = delta,
  f_t = f_t,
  dt_sub = dt
)

str(res)
```

This example performs a single infiltration–redistribution step consistent with the Struthers model equations, returning updated front positions, water contents, and drainage fluxes. Multiple steps can be iterated to represent continuous wet and dry cycles.

---

## Theoretical background

The Struthers model simplifies the Richards equation by assuming that soil moisture evolves as a sequence of sharp wetting fronts separating regions of nearly uniform water content. Infiltration creates new fronts, while redistribution gradually equalises moisture differences. The gravitational flux between adjacent layers is proportional to their hydraulic conductivity, and redistribution continues until the storage gradient vanishes or drainage occurs at the bottom boundary. Front merging ensures computational efficiency while maintaining mass balance.

The model can be applied to a wide range of soil textures by specifying \(K_s\), \(\theta_s\), \(\theta_r\), and \(\beta\) according to measured or literature-based hydraulic properties. Drainage can be parameterised by the total soil depth \(L\) and an effective lower boundary condition.

---

## Example scripts

The package includes several **fully reproducible examples** under the `examples/` directory (installed together with the package).  
These scripts illustrate the use of the *Struthers multiple–wetting–front* model under different boundary conditions and evaporation–transpiration regimes.

| Example file | Description |
|---------------|-------------|
| `drydown1.R` | Basic drydown simulation (no ET) with hourly redistribution and free drainage. |
| `drydown2.R` | Nonlinear drydown with daily evapotranspiration using the $\theta_\star$–$\theta_\mathrm{wp}$ stress formulation. |
| `example1.R` | Short infiltration–redistribution–drainage event with diagnostic plots. |
| `example2.R` | Full 40-day simulation combining rainfall bursts, ET, and depth-averaged moisture diagnostics. |

To **locate** these example files in your local installation, use `system.file()`:

```r
# List all example scripts available in the installed package
list.files(system.file("examples", package = "strutR"))

# Open one in the R editor
file.edit(system.file("examples", "drydown1.R", package = "strutR"))
```

You can also **run them directly** without copying them manually:

```r
source(system.file("examples", "drydown2.R", package = "strutR"))
```

Each example produces diagnostic plots and prints intermediate results,  
making it suitable both for research reproducibility and for teaching purposes.

---

## Citation

If you use this software, please cite **both** the original model and this R implementation:

Struthers, I., Hinz, C., & Sivapalan, M. (2006). *A multiple wetting front gravitational infiltration and redistribution model for water balance applications.* Water Resources Research, 42(6), W06406. https://doi.org/10.1029/2005WR004645

Martini, T. (2025). *strutR: Struthers Wetting-Front Model for Infiltration, Redistribution, and Drainage (R package).* Version 0.1.0. MIT License. https://doi.org/10.5281/zenodo.xxxxxx

---

## License and authorship

Copyright (c) 2025
Tommaso Martini  
Department of Regional and Urban Studies and Planning (DIST), University of Turin, Italy

Released under the MIT License. See [LICENSE](LICENSE) for details.

This software is an independent, open-source reimplementation of the Struthers model. It is **not affiliated** with the original authors or the American Geophysical Union (AGU), and it is distributed for research and educational purposes only.

---

## Related resources

- Original model publication: https://doi.org/10.1029/2005WR004645  
- GitHub repository: https://github.com/tommaso-martini/strutR  
- Example vignettes: `vignette("struthers-examples", package = "strutR")`  
