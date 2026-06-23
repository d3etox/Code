[🇫🇷 Version française](README.md)

# Uncertainty quantification & sensitivity analysis (R)

> Global sensitivity analysis and uncertainty propagation on the Borehole model, using Monte Carlo methods and Sobol indices.

**Context**: B. Iooss's sensitivity-analysis course (Polytech Nice-Sophia, 2026).

## What it does
On the reference **Borehole** model (water flow through a borehole, 8 input variables):
- Uncertainty propagation via **Monte Carlo** (convergence study, histograms)
- **Sobol indices**: variance decomposition to rank the influence of input parameters
- Scatterplots and graphical diagnostics

Skills: statistics, uncertainty quantification (UQ), R.

## Stack
**R** (RStudio).

## Files
- `borehole.R` — full analysis (author: Charles Fabre)
- `TP4-*.R`, `TP5sol-*.R` — course labs
- `*.png`, `Rplots.pdf` — results (MC convergence, Sobol indices…)

## Run
```r
source("borehole.R")
```
