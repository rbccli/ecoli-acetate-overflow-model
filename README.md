# Supplementary code for **NADH regulation shapes competing acetate flux responses in Escherichia coli : kinetic modelling of the tricarboxylic acid cycle and acetate production**.

The repository contains R/Stan code for Bayesian inference and MATLAB code for kinetic modelling, identifiability, sensitivity and regulation analyses.

## Code guide

All paths below are relative to `src/`.

| Purpose | Code |
| --- | --- |
| MATLAB model and parameters | `matlab/global_parameters_balanced.m`, `matlab/tca_reaction_fluxes_model/` |
| Structural identifiability | `matlab/identifiability analysis/tca_demo_genssi.m`, `matlab/identifiability analysis/tca_demo_genssi_run.m` |
| Profile likelihood for synthetic and experimental observations (Figures 10 and 13): shared algorithm, different inputs | `matlab/identifiability analysis/profile_likelihood.m` |
| Conservation analysis | `matlab/conserved_quantity_stability_analysis/conserved_quantity.m` |
| Bayesian calibration with experimental data | `r_stan/tcaferm_baye_real.R`, `r_stan/tcaferm_real.stan` |
| Experimental posterior predictive checks, parameter distributions and MCMC diagnostics | `r_stan/tcaferm_ppc_real.R`, `r_stan/tcaferm_prior_posterior_compare.R`, `r_stan/Diagnosis.R` |
| Total-order Sobol sensitivity analysis | `matlab/sensitivity analysis/sobol_sensitivity_total_order.m` |
| Regulation scenarios under redox-ratio and NAD(H)-pool perturbations | `matlab/redox compare/fourscenario.m`, `matlab/redox compare/fourscenario1.m` |
| Partial NADH inhibition and fermentation-response analysis | `matlab/redox compare/beta_inhibition.m`, `matlab/redox compare/beta_fermentation_two_perturbations.m` |
| Prior predictive simulation and prior/posterior predictive comparison | `r_stan/run_fast_stan_prior_predictive.R`, `r_stan/tca_ferm_prior_rng_fixed_param.stan`, `r_stan/tca_ferm_prior_posterior_predictive.stan` |
| Synthetic-data generation, Bayesian fitting and posterior predictive checks | `r_stan/synthetic/` |
| Shared figure formatting and export functions | `matlab/plot_helpers/`, `r_stan/recovered_plot_helpers.R` |

## Software dependencies

- **R:** CmdStan via `cmdstanr`, plus the R packages loaded by `library()` in the analysis scripts you run.
- **MATLAB:** Optimization, Parallel Computing, Statistics and Machine Learning, and Symbolic Math Toolboxes.
- **GenSSI:** required for structural identifiability. Add GenSSI and the relevant `src/matlab/` subdirectories to the MATLAB path.
