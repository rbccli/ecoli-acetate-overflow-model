# Supplementary code

This repository accompanies the manuscript **“NADH regulation shapes competing acetate flux responses in *Escherichia coli*: kinetic modelling of the tricarboxylic acid cycle and acetate production.”** It contains R/Stan code for Bayesian inference and MATLAB code for kinetic modelling, identifiability, sensitivity and regulation analyses.

## Code guide

All paths below are relative to `src/`.

| Analysis | Main files |
| --- | --- |
| Kinetic model and parameters | `matlab/global_parameters_balanced.m`; `matlab/tca_reaction_fluxes_model/` |
| Structural identifiability and profile likelihood | `matlab/identifiability analysis/` |
| Conservation and sensitivity analyses | `matlab/conserved_quantity_stability_analysis/`; `matlab/sensitivity analysis/` |
| Redox and NAD(H)-pool perturbations | `matlab/redox compare/` |
| Bayesian calibration and model checks using experimental observations | `r_stan/tcaferm_baye_real.R`; `r_stan/tcaferm_ppc_real.R`; `r_stan/tcaferm_prior_posterior_compare.R`; `r_stan/Diagnosis.R` |
| Prior predictive analyses | `r_stan/prior_posterior_predictive.R`; `r_stan/run_fast_stan_prior_predictive.R` |
| Synthetic-data analyses | `r_stan/synthetic/` |

## Requirements

- **R:** `cmdstanr`, CmdStan and a compatible C++ toolchain, plus the packages loaded by the R scripts.
- **MATLAB:** Optimization, Parallel Computing, Statistics and Machine Learning, and Symbolic Math Toolboxes.
- **GenSSI:** required for structural identifiability analysis.

## Getting started

Run the R analyses from src/r_stan/. Start with tcaferm_baye_real.R to fit the model to the experimental observations. The posterior predictive checks, parameter comparisons and diagnostics use the fitted model from that analysis.

For MATLAB analyses, add src/matlab/ and its subdirectories to the MATLAB path, then run the relevant script from its own directory. Some regulation scripts require posterior samples exported from the R fit and a reference-baseline file.

The scripts contain the observations and parameter settings used in the analyses; generated results and figures are created when the scripts are run.
