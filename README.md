# Supplementary code

This repository accompanies the manuscript:

**NADH regulation shapes competing acetate flux responses in *Escherichia coli*: kinetic modelling of the tricarboxylic acid cycle and acetate production**

It contains R/Stan code for Bayesian inference and MATLAB code for kinetic modelling, identifiability, sensitivity and regulation analyses.

## Code guide

All paths in this table are relative to `src/`.

| Purpose | Code |
| --- | --- |
| MATLAB model and baseline parameters | `matlab/global_parameters_balanced.m`, `matlab/tca_reaction_fluxes_model/` |
| Structural identifiability | `matlab/identifiability analysis/tca_demo_genssi.m`, `matlab/identifiability analysis/tca_demo_genssi_run.m` |
| Profile likelihood for synthetic and experimental observations (Figures 10 and 13; shared algorithm with separately configured inputs) | `matlab/identifiability analysis/profile_likelihood.m` |
| Conservation analysis | `matlab/conserved_quantity_stability_analysis/conserved_quantity.m` |
| Bayesian calibration with experimental observations | `r_stan/tcaferm_baye_real.R`, `r_stan/tcaferm_real.stan` |
| Experimental posterior predictive checks | `r_stan/tcaferm_ppc_real.R` |
| Prior/posterior parameter distributions and MCMC diagnostics | `r_stan/tcaferm_prior_posterior_compare.R`, `r_stan/Diagnosis.R` |
| Total-order Sobol sensitivity analysis | `matlab/sensitivity analysis/sobol_sensitivity_total_order.m` |
| Regulation scenarios under redox-ratio and NAD(H)-pool perturbations | `matlab/redox compare/fourscenario.m`, `matlab/redox compare/fourscenario1.m` |
| Partial NADH inhibition and fermentation-response analysis | `matlab/redox compare/beta_inhibition.m`, `matlab/redox compare/beta_fermentation_two_perturbations.m` |
| R-based prior predictive simulation and prior/posterior predictive comparison | `r_stan/prior_posterior_predictive.R` |
| Stan-based prior predictive simulation and prior/posterior predictive comparison | `r_stan/run_fast_stan_prior_predictive.R`, `r_stan/tca_ferm_prior_rng_fixed_param.stan`, `r_stan/tca_ferm_prior_posterior_predictive.stan` |
| Synthetic-data generation, Bayesian fitting and posterior predictive checks | `r_stan/synthetic/` |
| Shared figure formatting and export functions | `matlab/plot_helpers/`, `r_stan/recovered_plot_helpers.R` |

## Software requirements

### R and Stan

The analysis scripts use the following R packages:

`cmdstanr`, `jsonlite`, `dplyr`, `tidyr`, `ggplot2`, `posterior`, `patchwork`, `scales`, `bayesplot`, `deSolve`, `gtable` and `tibble`.

CmdStan and a compatible C++ toolchain are also required to compile and run the Stan models. Follow the [official CmdStanR installation instructions](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

The fitting scripts include an `install_cmdstan()` call. Once CmdStan is installed and configured, this installation step can be skipped.

### MATLAB

Depending on the analysis, the scripts require:

- Optimization Toolbox
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox

GenSSI is additionally required for structural identifiability analysis.

### Software versions

[TO COMPLETE: record the R, MATLAB, CmdStan, cmdstanr and GenSSI versions used to generate the manuscript results.]

## Input data and parameters

Experimental observations are embedded in `src/r_stan/tcaferm_baye_real.R`; the corresponding observed values are also included in the posterior predictive plotting script.

MATLAB baseline parameters and initial concentrations are defined in `src/matlab/global_parameters_balanced.m`.

The synthetic-data fitting script contains a fixed observation set. The data-generation script does not automatically replace these observations. To fit a newly generated dataset, update `y_obs`, `f_obs_ferm` and `f_obs_tca` in `tcaferm_baye_simulation.R` accordingly.

Posterior samples, fitted-model objects and reference-baseline files are generated during analysis and are not included in the repository.

## Running the R/Stan analyses

### Experimental-data calibration and post-processing

Start R from the repository root, then run:

```r
setwd("src/r_stan")

DEFAULT_FIG_PATH <- file.path(getwd(), "figures")
source("recovered_plot_helpers.R")

source("tcaferm_baye_real.R")
source("tcaferm_ppc_real.R")
source("tcaferm_prior_posterior_compare.R")
source("Diagnosis.R")
```

Run these scripts in the same R session: the post-processing scripts use the `tca_ferm_fit_real` object created by the fitting script.

The current fitting script performs both a full fit (`tca_ferm_fit_real`) and an additional shorter fit (`tca_ferm_fit_real_test`). The post-processing scripts listed above use the full fit.

The fitting script saves the full fit, posterior draws, summary statistics and CmdStan output files under:

```text
src/r_stan/results/stan_subunit_<timestamp>/
```

To reload a previously saved fit for posterior predictive checks or diagnostics:

```r
library(cmdstanr)

tca_ferm_fit_real <- readRDS(
  "results/stan_subunit_<timestamp>/stan_fit.rds"
)
```

Replace `<timestamp>` with the actual results directory. Retain the associated CmdStan output files, particularly when running generated-quantities analyses.

### Prior predictive analyses

After experimental-data calibration, run the selected prior predictive workflow:

```r
source("prior_posterior_predictive.R")
```

or:

```r
source("run_fast_stan_prior_predictive.R")
```

Both workflows require `data_list` and `tca_ferm_fit_real` in the R environment. Reloading `stan_fit.rds` alone does not recreate `data_list`.

The R-based workflow saves simulations, predictive summaries and figures under `results/prior_predictive_<timestamp>/`. The Stan-based workflow saves `stan_prior_predictive_draws_fast.rds` in the working directory and uses the shared plotting helpers for figure export.

### Synthetic-data analyses

Use a separate R session, starting from the repository root:

```r
setwd("src/r_stan/synthetic")

library(patchwork)

DEFAULT_FIG_PATH <- file.path(getwd(), "figures")
source("../recovered_plot_helpers.R")

source("tcaferm_baye_simulation.R")
source("tcaferm_ppc_simulation.R")
```

The posterior predictive script uses the `tcaferm_fit_simulation` object created by the fitting script.

To generate synthetic observations, run `tcaferm_simulation_estimation_scientific.R`. Transfer the selected generated observations into the fitting script before fitting a new synthetic dataset.

## Preparing posterior samples for MATLAB

The MATLAB regulation scripts read posterior parameters by column position. They require the 15 transformed model parameters in this order:

```text
PDH, ACO, AKGDH, SCS, SDH, PK, CS, IDH, PTA,
acetate, MDH, fumarase, ACKA, k, alpha
```

The general CSV export at the end of `tcaferm_baye_real.R` includes additional Stan variables. For MATLAB input, replace that export with the following ordered, theta-only export after fitting:

```r
theta_draws <- posterior::as_draws_matrix(
  tca_ferm_fit_real$draws("theta")
)

theta_columns <- sprintf("theta[%d]", seq_len(15))

write.csv(
  theta_draws[, theta_columns, drop = FALSE],
  file = "posterior1024.csv",
  row.names = FALSE
)
```

Place this file in the MATLAB analysis working directory, or update the relevant `readtable()` call to its location.

Most regulation scripts currently read `posterior1024.csv`; `beta_inhibition.m` reads `posterior0915.csv`. Update the filename to refer to the intended ordered posterior sample file before running the analysis.

## Running the MATLAB analyses

From the repository root, add the MATLAB source directories to the search path:

```matlab
addpath(genpath(fullfile(pwd, 'src', 'matlab')));
```

Then set the MATLAB Current Folder to the directory containing the selected analysis script and run that script. Relative input and output filenames are resolved from the working directory.

### Structural identifiability

Add GenSSI to the MATLAB path, then run:

```matlab
tca_demo_genssi_run
```

The runner loads the symbolic model defined in `tca_demo_genssi.m`.

### Profile likelihood

Run `profile_likelihood.m` separately for each observation set.

The checked-in script constructs `data.values` from concentrations in `global_parameters_balanced.m` and two flux observations. It does not provide an automatic switch between synthetic and experimental datasets.

For each dataset:

1. Set `data.values` to the corresponding observations, preserving the script's ordering of 13 concentrations followed by the two flux observations.
2. Check `data.obs_idx`, `y0`, `tspan`, `sigma_s` and `sigma_f` against the intended analysis.
3. Use distinct output filenames for each dataset to avoid overwriting results.

The script saves `profile_likelihood_results.mat` and exports the profile-likelihood figure.

### Sensitivity and regulation analyses

Run the selected script listed in the code guide after configuring its inputs.

For `fourscenario.m` and `fourscenario1.m`, first configure the PDH + CS reference scenario and set:

```matlab
compute_and_save_reference = true;
```

This creates the corresponding reference-baseline file. For subsequent scenario comparisons, set the flag to `false` so the scripts reuse the same posterior sample rows and paired reference baselines.

The reference files are:

- `posterior_reference_baselines_alpha.mat` for `fourscenario.m`
- `posterior_reference_baselines.mat` for `fourscenario1.m`

Scenario settings and model functions must be configured consistently; the scripts do not automatically generate all four scenarios in a single run.

## Outputs and reproducibility

R plotting scripts using `recovered_plot_helpers.R` save figures to `DEFAULT_FIG_PATH`. The MATLAB `export_figure()` helper saves figures to a `figures/` subdirectory of the working directory by default. Some MATLAB figure-export calls are commented out and must be enabled to save those figures.

Bayesian fitting, profile likelihood and global sensitivity analysis can require substantial computation. Sampling settings, simulation grids, solver tolerances and parallel settings are specified within the individual scripts.

When regenerating manuscript results, retain the input settings, random seeds, posterior sample files and software-version information associated with each analysis.

## License

[TO COMPLETE: specify the code license and add the corresponding LICENSE file to the repository.]
