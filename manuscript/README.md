# Sierra Nevada Bat Activity

This repository contains the data-processing, Bayesian modeling,
diagnostic, visualization, and manuscript workflow for an analysis of
seasonal and environmental drivers of bat acoustic activity in the
northern Sierra Nevada of California.

The study uses long-duration passive acoustic monitoring from eight
sites surveyed from April through October 2015. The analysis separates
broad seasonal phenology from nightly environmental associations and
evaluates whether responses to relative humidity differ among bat
assemblages characterized by low-, mid-, and high-frequency echolocation
calls.

## Scientific objectives

The current analysis addresses four related questions:

1.  How does bat activity vary over the active season?
2.  How are nightly temperature, relative humidity, and lunar phase
    associated with bat activity after accounting for seasonal
    phenology?
3.  Does humidity sensitivity differ among characteristic
    echolocation-frequency groups in the manner predicted by atmospheric
    attenuation theory?
4.  How much residual temporal dependence remains after accounting for
    site-specific seasonal phenology and environmental conditions?

A directed acyclic graph (DAG) is used to make the assumed relationships
among day of year, temperature, relative humidity, lunar phase, and bat
activity explicit and to guide model specification and interpretation.

## Analytical dataset

Detector effort is reconstructed before modeling so that nights with
detector failure are distinguished from sampled nights with zero
detected activity. The corrected analytical dataset contains 1,056
site-nights across eight sites.

Bat activity is summarized as the nightly number of 5-minute intervals
containing at least one classified bat call. Analyses are conducted for
the complete assemblage and for three characteristic
echolocation-frequency groups:

-   **Low frequency:** \<25 kHz
-   **Mid frequency:** 25--40 kHz
-   **High frequency:** \>40 kHz

The frequency groups are used to test predictions arising from the
frequency dependence of atmospheric attenuation.

## Modeling framework

The final analyses use Bayesian negative-binomial generalized additive
models. Seasonal phenology is represented flexibly and allowed to vary
among sites so that local seasonal trajectories are not forced into a
single regional pattern.

The manuscript distinguishes two complementary quantities:

-   **Total seasonal pattern:** activity across day of year under the
    environmental conditions that naturally occurred during the study.
-   **Conditional environmental model:** associations of temperature,
    relative humidity, lunar phase, and residual seasonal structure when
    the other modeled environmental variables are held constant.

Site-specific seasonal smooths are included because simpler seasonal
structures left substantial residual temporal dependence. Calendar-aware
residual autocorrelation diagnostics are used so that gaps caused by
detector failure or other missing sampling nights are not treated as
consecutive observations.

Model development emphasizes adequate representation of the ecological
and sampling processes relevant to the study questions rather than
adding complexity solely to reproduce every feature of the observed
data.

## Project structure

``` text
.
├── _quarto.yml
├── manuscript.qmd
├── sections/
│   ├── abstract.qmd
│   ├── introduction.qmd
│   ├── methods.qmd
│   ├── results.qmd
│   └── discussion.qmd
├── supplementary-material.qmd
├── references/
│   ├── references.bib
│   └── packages.bib
├── data/
├── derived/
├── figures/
├── tables/
├── scripts/
├── stan_programs/
└── _output/
```

The exact analysis-script names may evolve as the workflow is finalized.
Analysis scripts should generate derived data, posterior summaries,
tables, and figures; manuscript files should consume those outputs
rather than duplicate the analysis.

## Quarto manuscript workflow

The manuscript is organized as modular Quarto files. The master
`manuscript.qmd` should include the section files in manuscript order.
Section files are intended to be fragments of the complete manuscript
and therefore should not contain independent YAML headers when included
by the master document.

Project-level rendering defaults are stored in `_quarto.yml`.

From the project root, render the manuscript with:

``` bash
quarto render manuscript.qmd
```

The current project configuration writes rendered output to `_output/`.

The supplementary material can be rendered separately if it is
maintained as a standalone Quarto document:

``` bash
quarto render supplementary-material.qmd
```

## Analysis-to-manuscript workflow

The intended workflow is:

1.  Reconstruct and verify detector effort.
2.  Build the corrected aggregate and frequency-group nightly datasets.
3.  Fit candidate and final Bayesian GAMs.
4.  Evaluate sampling diagnostics, posterior predictive checks, and
    calendar-aware residual autocorrelation.
5.  Generate manuscript-ready posterior summaries, tables, and figures
    from the final model objects.
6.  Render the manuscript and supplementary material from those
    generated outputs.

Exploratory or model-development objects should not be treated as final
manuscript results. Numerical results, tables, and figures should be
regenerated from the accepted production fits.

## Current status

The aggregate model structure has been developed and evaluated,
including detector-effort correction, negative-binomial observation
models, site-specific seasonal phenology, posterior predictive
assessment, and residual temporal diagnostics.

Frequency-group models are being evaluated under the same framework.
Interpretations that depend on the frequency-group models---particularly
the form of humidity and temperature responses---should be considered
provisional until those production fits and diagnostics are complete.

The manuscript has working Quarto drafts for the Abstract, Introduction,
Methods, Results, Discussion, and Supplementary Material. Text inherited
from earlier model versions is being updated as final production-model
results become available.

## Reproducibility

Final model fits should be generated with multiple MCMC chains and
checked for convergence and sampling quality using R-hat, effective
sample sizes, divergent transitions, and related diagnostics. Posterior
predictive checks and residual temporal diagnostics are retained because
satisfactory MCMC sampling alone does not establish ecological model
adequacy.

Software citations and package versions should be generated from the
final computational environment and included with the manuscript
materials.

## Repository notes

Large fitted model objects, temporary outputs, rendered manuscript
files, and other reproducible intermediate products may be excluded from
version control through `.gitignore`. Source data, analysis code, model
programs, manuscript source files, and lightweight derived summaries
needed to reproduce the reported results should be retained according to
project data-sharing requirements.

## Citation

A formal citation and data/code availability statement will be added
when the manuscript and repository are prepared for submission or
archival release.
