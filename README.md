# Solar Energy Risks Replication Package

This replication package accompanies the paper *Solar Energy Risks: Stochastic
Radiation Modelling and Optimal Hedging Strategies*.

The citable archive is available on Zenodo:
<https://doi.org/10.5281/zenodo.18065941>.

## Archive and repository policy

The replication material is intended to be distributed in two forms.

1. A fixed Zenodo snapshot, used for citation and long-term preservation.
2. A selective public Git mirror containing the README, `outputs.RData`, the
   small public input/example data files, and the scripts/functions needed to
   reproduce the empirical pipeline.

Private submission files and revision-management material, such as the main
manuscript, appendix, `revision/`, and rendered site sources, are maintained
outside the selective public mirror. The full archival snapshot on Zenodo may
contain additional rendered documentation and environment files.

## Folder structure

- `scripts/`: R scripts used to reproduce the empirical pipeline.
- `scripts/data/`: model estimation, simulations, pricing inputs, moments,
  diagnostics, and forecasts.
- `scripts/tables/`: table-generation scripts for the manuscript, appendix,
  and Supplementary Material.
- `scripts/figs/`: figure-generation scripts.
- `scripts/functions/`: local helper functions for the electricity model,
  radiation model, CTMC calculations, simulations, pricing, figures, and file
  output.
- `scripts/tests/testthat/`: unit tests for CTMC and radiation-model helper
  functions.
- `data/`: input data and generated model, simulation, pricing, and diagnostic
  objects. In the selective public mirror, only `data/GME_daily.csv` and the
  curated `data/exm/` example objects are included; remaining generated objects
  should be restored from the Zenodo snapshot before rerunning the full
  pipeline.
- `figs/`: generated figures used in the manuscript and appendix.
- `outputs.RData`: central registry of output paths, settings, and generated
  table objects.
- `environment/`: reproducible R environment files, including `renv.lock` and
  the full package-version snapshot CSV, included in the full archival snapshot
  when available.
- `main.qmd`, `appendix.qmd`, `SM1-supplementary-material.qmd`: manuscript,
  appendix, and Supplementary Material sources, kept outside the selective
  public mirror unless explicitly distributed with the full archive.
- `site_/`: source and rendered files for the companion website, kept outside
  the selective public mirror.

## Software requirements

The scripts were prepared and tested with R 4.4.1. They should be run from the
root of the replication folder.

The R dependency audit was performed on the active replication sources
(`scripts/`, `scripts/tests/`, the manuscript/appendix/Supplementary Material
sources, and the companion-site sources when present). Deprecated folders,
revision-only scripts, and generated HTML folders are not part of the
replication dependency set.

Core R dependencies are:

`solarr`, `readxl`, `Rcpp`, `tidyverse`, `purrr`, `dplyr`, `readr`, `tibble`,
`tidyr`, `stringr`, `ggplot2`, `lubridate`, `broom`, `R6`, `mixtools`,
`numDeriv`, `expm`, `cubature`, `mvtnorm`, `quadprog`, `quantreg`, `cli`,
`crayon`, `gridExtra`, `gridtext`, and `latex2exp`.

For tables, tests, animations, spell checks, and rendered documents, the
additional requirements are:

`knitr`, `kableExtra`, `testthat`, `gganimate`, `spelling`, `renv`, and Quarto.

The package `solarr` is the project package used by the replication scripts. It
is available from GitHub at <https://github.com/beniamino98/solarr>. In the
current lockfile it is recorded as version `1.0.1`. If `renv::restore()` cannot
retrieve it automatically on a new machine, install it from GitHub first and
then run the restore command below:

```r
install.packages("remotes")
remotes::install_github("beniamino98/solarr")
```

To reproduce the package versions used for the full snapshot, use the included
lockfile under `environment/renv.lock` when it is available:

```r
install.packages("renv")
renv::restore(lockfile = "environment/renv.lock", prompt = FALSE)
```

For quick inspection, `environment/r-package-snapshot.csv`, when available,
contains the package names, versions, source type, and repository recorded from
the current R library.

The folder `scripts/functions/C/` contains C source and compiled shared objects
used by some numerical routines. Recompilation requires a working C toolchain
compatible with the local R installation.

## Replication pipeline

The pipeline is organized in four top-level stages. Use `RScript` if that alias
is available on the system; otherwise replace it with `Rscript`. The orchestration
scripts also respect the optional environment variable `RSCRIPT_BIN`.

```bash
RScript scripts/s0-outputs.R
RScript scripts/s1-data.R "TRUE" "5000" "2022" "1,2,3,5,10,15,30"
RScript scripts/s2-tables.R "TRUE"
RScript scripts/s3-figures.R "TRUE" "2022" "2014" "FALSE"
```

The stages are:

1. `scripts/s0-outputs.R`
   initializes `outputs.RData` and the expected output directories.

2. `scripts/s1-data.R`
   runs the complete data-generation pipeline. It estimates the radiation
   models under the historical probability measure, estimates electricity-price
   models and pricing inputs, estimates residual correlations, simulates joint
   radiation-electricity paths, computes contract moments and mean-variance
   objects, and generates diagnostic and forecast data for the Supplementary
   Material.

3. `scripts/s2-tables.R`
   reads the generated objects and stores manuscript, appendix, and
   Supplementary Material tables in `outputs.RData`.

4. `scripts/s3-figures.R`
   reads the generated objects and writes the figures under `figs/`.

The data stage calls the scripts in `scripts/data/` in this order:

```text
s0-models-radiation-P-discrete-place.R
s1-models-radiation-P-IID-place.R
s2a-models-radiation-P-DTMC-place.R
s2b-models-radiation-P-CTMC-place.R
s3-models-electricity-P.R
s4-models-electricity-Q.R
s5-models-rho-place.R
s6-simulate-Rt-Et-place.R
s7-solarOptions-moments-place.R
s8-solarOptions-mv-place.R
s9-models-radiation-P-moments-place.R
s10-bounds-2gm.R
SM1-models-radiation-P-diagnostic-place.R
SM2-models-radiation-P-forecasts-place.R
```

The default empirical settings are the three locations `Bologna`, `Roma`, and
`Palermo`, `5000` Monte Carlo scenarios, reference year `2022` for the
two-Gaussian bounds, and horizons `1, 2, 3, 5, 10, 15, 30`.

## Companion website

The companion website in `site_/`, when available in the full archive,
documents the same pipeline at script level and provides the compiled
Supplementary Material. It is not required for rerunning the statistical
pipeline, and it is not part of the selective public mirror.
