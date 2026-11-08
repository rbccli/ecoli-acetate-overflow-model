# ============================================================
# Prior predictive simulation aligned with the real-data Stan model
# Requires data_list and the saved tca_ferm_fit_real; does not run MCMC.
# options(tca.prior.n_draws = 20L) can be used for a short smoke run.
# ============================================================

library(deSolve)
library(posterior)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(20260816)

# ============================================================
# 1. Validate the actual Stan data and load its priors/transforms
# ============================================================

if (!exists("data_list", inherits = TRUE)) {
  stop("Load the real-data data_list before running this script.", call. = FALSE)
}
if (!exists("tca_ferm_fit_real", inherits = TRUE)) {
  stop("Load the saved stan_fit.rds as tca_ferm_fit_real first.", call. = FALSE)
}

# When sourced by full path, resolve the model beside this script. When run
# interactively in RStudio, use the working directory or an explicit option.
source_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
script_dir <- if (length(source_files)) dirname(tail(source_files, 1L)) else getwd()
stan_model_file <- getOption(
  "tca.prior.stan_file",
  file.path(script_dir, "tcaferm_real.stan")
)
if (!file.exists(stan_model_file)) {
  stop("Cannot find tcaferm_real.stan. Set options(tca.prior.stan_file = '/full/path/tcaferm_real.stan').",
       call. = FALSE)
}
stan_model_file <- normalizePath(stan_model_file, mustWork = TRUE)
stan_lines <- sub("//.*$", "", readLines(stan_model_file, warn = FALSE))

prior_matches <- regmatches(
  stan_lines,
  regexec("^\\s*(\\w+_raw)\\s*~\\s*normal\\s*\\((.*)\\)\\s*;", stan_lines, perl = TRUE)
)
prior_matches <- Filter(function(x) length(x) == 3L, prior_matches)
prior_spec <- do.call(rbind, lapply(prior_matches, function(x) {
  values <- eval(str2lang(paste0("c(", x[3], ")")), envir = baseenv())
  stopifnot(length(values) == 2L, all(is.finite(values)), values[2] > 0)
  data.frame(variable = x[2], mean = values[1], sd = values[2])
}))
expected_raw <- c("pdh_raw", "aco_raw", "akd_raw", "scs_raw", "sd_raw",
                  "pk_raw", "cs_raw", "id_raw", "pta_raw", "acetate_raw",
                  "md_raw", "f_raw", "ack_raw", "k_raw", "alpha_raw",
                  "sigma_c_raw", "sigma_f_raw")
if (is.null(prior_spec) || anyDuplicated(prior_spec$variable) ||
    !setequal(prior_spec$variable, expected_raw)) {
  stop("The Stan raw-parameter priors differ from the supported real-data model.", call. = FALSE)
}
prior_spec <- prior_spec[match(expected_raw, prior_spec$variable), , drop = FALSE]
rownames(prior_spec) <- NULL

# Read the exact theta ordering and fixed acetate from transformed parameters.
theta_matches <- regmatches(
  stan_lines,
  regexec("^\\s*theta\\[([0-9]+)\\]\\s*=\\s*([^;]+);", stan_lines, perl = TRUE)
)
theta_matches <- Filter(function(x) length(x) == 3L, theta_matches)
theta_indices <- as.integer(vapply(theta_matches, `[[`, character(1), 2L))
if (!identical(theta_indices, seq_len(15L))) {
  stop("Unexpected theta order in tcaferm_real.stan.", call. = FALSE)
}
theta_expressions <- lapply(theta_matches, function(x) str2lang(x[3]))
theta_parameter_names <- c("enzyme_pdh", "enzyme_aco", "enzyme_akd", "enzyme_scs",
                           "enzyme_sd", "enzyme_pk", "enzyme_cs", "enzyme_id",
                           "enzyme_pta", "acetate", "enzyme_md", "enzyme_fumarase",
                           "enzyme_ack", "k", "alpha")
raw_draw_to_theta <- function(draw) {
  values <- as.list(draw)
  values$pow <- function(x, y) x^y
  vapply(theta_expressions, function(expr) {
    as.numeric(eval(expr, envir = values, enclos = baseenv()))
  }, numeric(1))
}

# Map the packed parameters actually supplied to Stan, rather than relying on
# an old params object or the synthetic-data script's hard-coded constants.
parameter_names <- list(
  p_pdh = c("km_co2_pdh", "km_pyruvate_pdh", "km_coa_pdh", "km_nad_pdh",
            "km_acetyl_coa_pdh", "km_nadh_pdh", "kf_pdh", "ki_nadh_pdh", "kr_pdh"),
  p_cs = c("km_oxaloacetate_cs", "km_acetyl_coa_cs", "km_citrate_cs", "km_coa_sh_cs",
           "kf_cs", "kr_cs", "ki_nadh_cs"),
  p_aco = c("km_citrate_aco", "km_cis_aconitate_aco1", "km_cis_aconitate_aco2",
            "km_isocitrate_aco", "kf1_aco", "kr1_aco", "kf2_aco", "kr2_aco"),
  p_id = c("km_co2_id", "km_isocitrate_id", "km_nadh_id", "km_nad_id",
           "km_aketo_id", "kf_id", "kr_id"),
  p_akd = c("km_co2_akd", "km_aketo_akd", "km_nad_akd", "km_nadh_akd",
            "km_coa_sh_akd", "km_succinyl_coa_akd", "kf_akd", "kr_akd"),
  p_scs = c("km_succinyl_coa_scs", "km_adp_scs", "km_succinate_scs", "km_coa_sh_scs",
            "km_atp_scs", "kf_scs", "kr_scs", "km_pi_scs"),
  p_sd = c("km_succinate_sd", "km_fad_sd", "km_fumarate_sd", "km_fadh2_sd", "kf_sd", "kr_sd"),
  p_fum = c("km_fumarate_fumarase", "km_malate_fumarase", "kf_fumarase", "kr_fumarase"),
  p_md = c("km_nad_md", "km_oxaloacetate_md", "km_malate_md", "km_nadh_md", "kf_md", "kr_md"),
  p_pk = c("km_pep_pk", "km_atp_pk", "km_adp_pk", "kf_pk", "km_pyruvate_pk", "ki_atp_pk", "kr_pk"),
  p_pta = c("km_acetyl_coa_pta", "km_pi_pta", "km_acetyl_p_pta", "km_coa_pta", "kf_pta", "kr_pta"),
  p_ack = c("km_acetyl_p_ack", "km_adp_ack", "km_acetate_ack", "km_atp_ack", "kf_ack", "kr_ack")
)
required_data <- c("NY", "N_t", "y0", "t_obs", "y_obs", "C", "f_obs_tca", "f_obs_ferm",
                   names(parameter_names))
if (!all(required_data %in% names(data_list))) {
  stop("data_list is missing: ", paste(setdiff(required_data, names(data_list)), collapse = ", "),
       call. = FALSE)
}
stopifnot(data_list$NY == 13L, data_list$N_t >= 1L,
          length(data_list$y0) == 13L, all(is.finite(data_list$y0)),
          length(data_list$C) == 11L, all(is.finite(data_list$C)),
          identical(dim(data_list$y_obs), c(as.integer(data_list$N_t), 13L)),
          length(data_list$t_obs) == data_list$N_t,
          all(is.finite(data_list$t_obs)), all(diff(c(0, data_list$t_obs)) > 0))
y0_prior <- as.numeric(data_list$y0)
times_prior <- c(0, as.numeric(data_list$t_obs))
base_params <- list(C = as.numeric(data_list$C))
for (block in names(parameter_names)) {
  values <- as.numeric(data_list[[block]])
  if (length(values) != length(parameter_names[[block]]) || any(!is.finite(values))) {
    stop("Invalid Stan parameter block: ", block, call. = FALSE)
  }
  base_params[parameter_names[[block]]] <- as.list(values)
}

# Match Stan's ode_bdf_tol budget/tolerances; VODE mf=22 is implicit BDF.
prior_rtol <- 1e-3
prior_atol <- 1e-5
prior_maxsteps <- 1000000L

# ============================================================
# 2. Draw parameters from the current Stan priors
# ============================================================
n_prior <- as.integer(getOption("tca.prior.n_draws", 5000L))
stopifnot(length(n_prior) == 1L, is.finite(n_prior), n_prior >= 1L)
prior_draws <- as.data.frame(setNames(lapply(seq_len(nrow(prior_spec)), function(j) {
  rnorm(n_prior, prior_spec$mean[j], prior_spec$sd[j])
}), prior_spec$variable))
cat("Stan prior source:", stan_model_file, "\n")
print(prior_spec, row.names = FALSE)

# ============================================================
# 3. Kinetics translated from this real-data Stan model
# ATP inhibition is disabled there; NADH inhibition is retained.
# ============================================================

atp_inhibition <- function(atp_con, ki_atp) {
  return(1.0)
}

nadh_inhibition <- function(nadh_con, ki_nadh) {
  return(1.0 / (1.0 + nadh_con / ki_nadh))
}

compute_v_aco <- function(citrate_con, cis_aconitate_con, isocitrate_con, km_cit, km_cis1, km_cis2, km_iso, enzyme_aco, kf1, kr1, kf2, kr2) {
  cit_tilde <- citrate_con / km_cit
  cis_tilde1 <- cis_aconitate_con / km_cis1
  cis_tilde2 <- cis_aconitate_con / km_cis2
  iso_tilde <- isocitrate_con / km_iso
  denom1 <- 1.0 + cit_tilde + cis_tilde1
  denom2 <- 1.0 + cis_tilde2 + iso_tilde
  v <- numeric(2)
  v[1] <- (kf1 * cit_tilde - kr1 * cis_tilde1) / denom1 * enzyme_aco
  v[2] <- (kf2 * cis_tilde2 - kr2 * iso_tilde) / denom2 * enzyme_aco
  return(v)
}

compute_v_akd <- function(a_keto_con, nad_con, coa_sh_con, succinyl_coa_con, nadh_con, co2_con, km_co2, km_ak, km_nad, km_nadh, km_coa, km_succinyl_coa, enzyme_akd, kf_akd, kr_akd) {
  ak_tilde <- a_keto_con / km_ak
  nad_tilde <- nad_con / km_nad
  coa_tilde <- coa_sh_con / km_coa
  succinyl_coa_tilde<- succinyl_coa_con / km_succinyl_coa
  nadh_tilde <- nadh_con / km_nadh
  co2_tilde <- co2_con / km_co2
  denom <- (1.0 + nad_tilde) * (1.0 + ak_tilde) * (1.0 + coa_tilde) + (1.0 + succinyl_coa_tilde) * (1.0 + nadh_tilde) * (1.0 + co2_tilde) - 1.0
  fwd <- kf_akd * ak_tilde * coa_tilde * nad_tilde
  rev <- kr_akd * succinyl_coa_tilde * nadh_tilde * co2_tilde
  v <- (fwd - rev) / denom * enzyme_akd
  return(v)
}

compute_v_cs <- function(oxaloacetate_con, acetyl_coa_con, citrate_con, coa_sh_con, nadh_con, km_ox, km_acc, km_cit, km_coa_sh, enzyme_cs, kf_cs, kr_cs, ki_nadh) {
  ox_tilde <- oxaloacetate_con / km_ox
  acc_tilde <- acetyl_coa_con / km_acc
  cit_tilde <- citrate_con / km_cit
  coa_tilde <- coa_sh_con / km_coa_sh
  denom <- (1.0 + ox_tilde) * (1.0 + acc_tilde) + (1.0 + cit_tilde) * (1.0 + coa_tilde) - 1.0
  fwd <- kf_cs * ox_tilde * acc_tilde
  rev <- kr_cs * cit_tilde * coa_tilde
  v <- (fwd - rev) / denom * enzyme_cs
  v <- v * (nadh_inhibition(nadh_con, ki_nadh))
  return(v)
}

compute_v_fumarase <- function(fumarate_con, malate_con, km_fum, km_mal, enzyme_fum, kf_fum, kr_fum) {
  fum_tilde <- fumarate_con / km_fum
  mal_tilde <- malate_con / km_mal
  denom <- 1.0 + fum_tilde + mal_tilde
  v <- (kf_fum * fum_tilde - kr_fum * mal_tilde) / denom * enzyme_fum
  return(v)
}

compute_v_id <- function(isocitrate_con, nad_con, nadh_con, aketo_con, co2_con, km_co2, km_iso, km_nadh, km_nad, km_aketo, enzyme_id, kf_id, kr_id) {
  iso_tilde <- isocitrate_con / km_iso
  nad_tilde <- nad_con / km_nad
  nadh_tilde <- nadh_con / km_nadh
  aketo_tilde <- aketo_con / km_aketo
  co2_tilde <- co2_con / km_co2
  denom <- (1.0 + iso_tilde) * (1.0 + nad_tilde) + (1.0 + nadh_tilde) * (1.0 + aketo_tilde) * (1.0 + co2_tilde) - 1.0
  fwd <- kf_id * iso_tilde * nad_tilde
  rev <- kr_id * nadh_tilde * aketo_tilde * co2_tilde
  v <- (fwd - rev) / denom * enzyme_id
  return(v)
}

compute_v_md <- function(malate_con, nad_con, oaa_con, nadh_con, km_mal, km_nad, km_oaa, km_nadh, enzyme_md, kf_md, kr_md) {
  mal_tilde <- malate_con / km_mal
  nad_tilde <- nad_con / km_nad
  oaa_tilde <- oaa_con / km_oaa
  nadh_tilde <- nadh_con / km_nadh
  denom <- (1.0 + mal_tilde) * (1.0 + nad_tilde) + (1.0 + oaa_tilde) * (1.0 + nadh_tilde) - 1.0
  v <- (kf_md * mal_tilde * nad_tilde - kr_md * oaa_tilde * nadh_tilde) / denom * enzyme_md
  return(v)
}

compute_v_pdh <- function(pyruvate_con, coa_con, nad_con, acetyl_coa_con, nadh_con, co2_con, km_co2, km_pyruvate, km_coa, km_nad, km_acetyl_coa, km_nadh, enzyme_pdh, kf_pdh, kr_pdh, ki_nadh) {
  pyruvate_tilde <- pyruvate_con / km_pyruvate
  coa_tilde <- coa_con / km_coa
  nad_tilde <- nad_con / km_nad
  acetyl_coa_tilde <- acetyl_coa_con / km_acetyl_coa
  nadh_tilde <- nadh_con / km_nadh
  co2_tilde <- co2_con / km_co2
  denom <- (1.0 + pyruvate_tilde) * (1.0 + coa_tilde) * (1.0 + nad_tilde) + (1.0 + acetyl_coa_tilde) * (1.0 + nadh_tilde) * (1.0 + co2_tilde) - 1.0
  fwd <- kf_pdh * pyruvate_tilde * coa_tilde * nad_tilde
  rev <- kr_pdh * acetyl_coa_tilde * nadh_tilde * co2_tilde
  v <- (fwd - rev) / denom * enzyme_pdh
  v <- v * (nadh_inhibition(nadh_con, ki_nadh))
  return(v)
}

compute_v_pk <- function(km_pyruvate, km_atp, km_pep, km_adp, kf_pk, kr_pk, enzyme_pk, pyruvate_con, atp_con, pep_con, adp_con, ki_atp) {
  pyruvate_tilde <- pyruvate_con / km_pyruvate
  pep_tilde <- pep_con / km_pep
  atp_tilde <- atp_con / km_atp
  adp_tilde <- adp_con / km_adp
  denom <- (1.0 + pyruvate_tilde) * (1.0 + atp_tilde) + (1.0 + pep_tilde) * (1.0 + adp_tilde) - 1.0
  v <- (kf_pk * pep_tilde * adp_tilde - kr_pk * pyruvate_tilde * atp_tilde) / denom * enzyme_pk
  return(v)
}

compute_v_scs <- function(succinyl_coa_con, adp_con, succinate_con, coa_sh_con, nadh_con, atp_con, pi_con, km_succinyl, km_adp, km_succinate, km_coa, km_atp, km_pi, enzyme_scs, kf_scs, kr_scs) {
  sucCoA_tilde <- succinyl_coa_con / km_succinyl
  adp_tilde <- adp_con / km_adp
  suc_tilde <- succinate_con / km_succinate
  coa_tilde <- coa_sh_con / km_coa
  atp_tilde <- atp_con / km_atp
  pi_tilde <- pi_con / km_pi
  denom <- (1.0 + sucCoA_tilde) * (1.0 + adp_tilde) * (1.0 + pi_tilde) + (1.0 + suc_tilde) * (1.0 + coa_tilde) * (1.0 + atp_tilde) - 1.0
  v <- (kf_scs * sucCoA_tilde * adp_tilde * pi_tilde - kr_scs * suc_tilde * coa_tilde * atp_tilde) / denom * enzyme_scs
  return(v)
}

compute_v_sd <- function(succinate_con, uq_con, fumarate_con, uqh2_con, km_suc, km_uq, km_fum, km_uqh2, enzyme_sd, kf_sd, kr_sd) {
  suc_tilde <- succinate_con / km_suc
  uq_tilde <- uq_con / km_uq
  fum_tilde <- fumarate_con / km_fum
  uqh2_tilde <- uqh2_con / km_uqh2
  denom <- (1.0 + suc_tilde) * (1.0 + uq_tilde) + (1.0 + fum_tilde) * (1.0 + uqh2_tilde) - 1.0
  v <- (kf_sd * suc_tilde * uq_tilde - kr_sd * fum_tilde * uqh2_tilde) / denom * enzyme_sd
  return(v)
}

compute_v_pta <- function(acetyl_coa_con, pi_con, acetyl_p_con, coa_sh_con, km_acetyl_coa_pta, km_pi_pta, km_acetyl_p_pta, km_coa_pta, kf_pta, kr_pta, enzyme_pta) {
  acetyl_coa_tilde <- acetyl_coa_con / km_acetyl_coa_pta
  pi_tilde <- pi_con / km_pi_pta
  acetyl_p_tilde <- acetyl_p_con / km_acetyl_p_pta
  coa_tilde <- coa_sh_con / km_coa_pta
  denom <- (1 + acetyl_coa_tilde) * (1 + pi_tilde) + (1 + acetyl_p_tilde) * (1 + coa_tilde) - 1
  v <- (kf_pta * acetyl_coa_tilde * pi_tilde - kr_pta * acetyl_p_tilde * coa_tilde)/ denom * enzyme_pta
  return(v)
}

compute_v_ack <- function(acetyl_p_con, adp_con, acetate_con, atp_con, km_acetyl_p_ack, km_adp_ack, km_acetate_ack, km_atp_ack, kf_ack, kr_ack, enzyme_ack) {
  acetyl_p_tilde <- acetyl_p_con / km_acetyl_p_ack
  adp_tilde <- adp_con / km_adp_ack
  atp_tilde <- atp_con / km_atp_ack
  acetate_tilde <- acetate_con / km_acetate_ack
  denom <- (1 + acetyl_p_tilde) * (1 + adp_tilde) + (1 + acetate_tilde) * (1 + atp_tilde) - 1
  v <- (kf_ack * acetyl_p_tilde * adp_tilde - kr_ack * acetate_tilde * atp_tilde) / denom * enzyme_ack
  return(v)
}

compute_v_redox <- function(nadh_con, nad_con, k, alpha) {
  v <- k * alpha * nad_con - k * nadh_con
  return(v)
}

compute_all_fluxes_R <- function(
    y,
    parms,
    C = parms$C,
    acetate = parms$acetate
) {
  
  pyruvate       <- y[1]
  oxaloacetate   <- y[2]
  acetyl_coa     <- y[3]
  citrate        <- y[4]
  cis_aconitate  <- y[5]
  isocitrate     <- y[6]
  aketo          <- y[7]
  succinyl_coa   <- y[8]
  succinate      <- y[9]
  fumarate       <- y[10]
  acetyl_p       <- y[11]
  nad            <- y[12]
  nadh           <- y[13]
  
  atp   <- C[3]
  adp   <- C[4]
  fad   <- C[5]
  fadh2 <- C[6]
  co2   <- C[7]
  pep   <- C[8]
  pi    <- C[11]
  
  coa_sh <- C[9] - acetyl_coa - succinyl_coa
  
  malate <- C[10] - (
    oxaloacetate +
      citrate +
      cis_aconitate +
      isocitrate +
      aketo +
      succinyl_coa +
      succinate +
      fumarate
  )
  
  with(parms, {
    
    v_pdh <- compute_v_pdh(
      pyruvate, coa_sh, nad,
      acetyl_coa, nadh, co2,
      km_co2_pdh,
      km_pyruvate_pdh,
      km_coa_pdh,
      km_nad_pdh,
      km_acetyl_coa_pdh,
      km_nadh_pdh,
      enzyme_pdh,
      kf_pdh,
      kr_pdh,
      ki_nadh_pdh
    )
    
    v_cs <- compute_v_cs(
      oxaloacetate,
      acetyl_coa,
      citrate,
      coa_sh,
      nadh,
      km_oxaloacetate_cs,
      km_acetyl_coa_cs,
      km_citrate_cs,
      km_coa_sh_cs,
      enzyme_cs,
      kf_cs,
      kr_cs,
      ki_nadh_cs
    )
    
    v_aco <- compute_v_aco(
      citrate,
      cis_aconitate,
      isocitrate,
      km_citrate_aco,
      km_cis_aconitate_aco1,
      km_cis_aconitate_aco2,
      km_isocitrate_aco,
      enzyme_aco,
      kf1_aco,
      kr1_aco,
      kf2_aco,
      kr2_aco
    )
    
    v_id <- compute_v_id(
      isocitrate,
      nad,
      nadh,
      aketo,
      co2,
      km_co2_id,
      km_isocitrate_id,
      km_nadh_id,
      km_nad_id,
      km_aketo_id,
      enzyme_id,
      kf_id,
      kr_id
    )
    
    v_akd <- compute_v_akd(
      aketo,
      nad,
      coa_sh,
      succinyl_coa,
      nadh,
      co2,
      km_co2_akd,
      km_aketo_akd,
      km_nad_akd,
      km_nadh_akd,
      km_coa_sh_akd,
      km_succinyl_coa_akd,
      enzyme_akd,
      kf_akd,
      kr_akd
    )
    
    v_scs <- compute_v_scs(
      succinyl_coa,
      adp,
      succinate,
      coa_sh,
      nadh,
      atp,
      pi,
      km_succinyl_coa_scs,
      km_adp_scs,
      km_succinate_scs,
      km_coa_sh_scs,
      km_atp_scs,
      km_pi_scs,
      enzyme_scs,
      kf_scs,
      kr_scs
    )
    
    v_sd <- compute_v_sd(
      succinate,
      fad,
      fumarate,
      fadh2,
      km_succinate_sd,
      km_fad_sd,
      km_fumarate_sd,
      km_fadh2_sd,
      enzyme_sd,
      kf_sd,
      kr_sd
    )
    
    v_fum <- compute_v_fumarase(
      fumarate,
      malate,
      km_fumarate_fumarase,
      km_malate_fumarase,
      enzyme_fumarase,
      kf_fumarase,
      kr_fumarase
    )
    
    v_md <- compute_v_md(
      malate,
      nad,
      oxaloacetate,
      nadh,
      km_malate_md,
      km_nad_md,
      km_oxaloacetate_md,
      km_nadh_md,
      enzyme_md,
      kf_md,
      kr_md
    )
    
    v_pk <- compute_v_pk(
      km_pyruvate_pk,
      km_atp_pk,
      km_pep_pk,
      km_adp_pk,
      kf_pk,
      kr_pk,
      enzyme_pk,
      pyruvate,
      atp,
      pep,
      adp,
      ki_atp_pk
    )
    
    v_pta <- compute_v_pta(
      acetyl_coa,
      pi,
      acetyl_p,
      coa_sh,
      km_acetyl_coa_pta,
      km_pi_pta,
      km_acetyl_p_pta,
      km_coa_pta,
      kf_pta,
      kr_pta,
      enzyme_pta
    )
    
    v_ack <- compute_v_ack(
      acetyl_p,
      adp,
      acetate,
      atp,
      km_acetyl_p_ack,
      km_adp_ack,
      km_acetate_ack,
      km_atp_ack,
      kf_ack,
      kr_ack,
      enzyme_ack
    )
    
    v_redox <- compute_v_redox(
      nadh,
      nad,
      k,
      alpha
    )
    
    c(
      PK = v_pk,
      PDH = v_pdh,
      CS = v_cs,
      ACO1 = v_aco[1],
      ACO2 = v_aco[2],
      IDH = v_id,
      AKGDH = v_akd,
      SCS = v_scs,
      SDH = v_sd,
      FUM = v_fum,
      MDH = v_md,
      PTA = v_pta,
      ACKA = v_ack,
      redox = v_redox
    )
  })
}


# One shared flux implementation feeds both the ODE and final flux output.
# State order: PYR, OAA, AcCoA, CIT, CIS, ISO, AKG, SucCoA, SUC, FUM,
# acetyl-P, NAD+, NADH. This is rhs() in tcaferm_real.stan.
ode_func <- function(t, y, parms) {
  v <- unname(compute_all_fluxes_R(y, parms))
  redox_balance <- v[2] + v[6] + v[7] + v[11] + v[14]
  list(c(
    v[1] - v[2],
    v[11] - v[3],
    v[2] - v[3] - v[12],
    v[3] - v[4],
    v[4] - v[5],
    v[5] - v[6],
    v[6] - v[7],
    v[7] - v[8],
    v[8] - v[9],
    v[9] - v[10],
    v[12] - v[13],
    -redox_balance,
    redox_balance
  ))
}

# ============================================================
# 4. Run one prior predictive simulation
# ============================================================

run_one_prior_draw <- function(i) {
  
  draw_i <- prior_draws[i, ]
  
  parms_i <- base_params
  theta_i <- raw_draw_to_theta(draw_i)
  parms_i[theta_parameter_names] <- as.list(theta_i)
  parms_i$theta <- theta_i

  sigma_c_i <- 10^draw_i$sigma_c_raw
  sigma_f_i <- 10^draw_i$sigma_f_raw
  
  ode_attempt <- tryCatch(
    {
      list(
        ok = TRUE,
        value = deSolve::ode(
          y = y0_prior,
          times = times_prior,
          func = ode_func,
          parms = parms_i,
          method = "vode",
          rtol = prior_rtol,
          atol = prior_atol,
          maxsteps = prior_maxsteps,
          mf = 22
        )
      )
    },
    error = function(e) {
      list(
        ok = FALSE,
        message = conditionMessage(e)
      )
    }
  )
  
  if (!isTRUE(ode_attempt$ok)) {
    return(
      list(
        success = FALSE,
        reason = paste0(
          "ODE failure: ",
          ode_attempt$message
        )
      )
    )
  }
  
  ode_output <- ode_attempt$value

  solver_status <- attr(ode_output, "istate")
  target_time <- tail(times_prior, 1L)
  reached_time <- tail(ode_output[, 1], 1L)
  if (length(solver_status) == 0L || solver_status[1] != 2L ||
      nrow(ode_output) != length(times_prior) ||
      !isTRUE(all.equal(as.numeric(reached_time), target_time, tolerance = 1e-10))) {
    return(list(success = FALSE, reason = paste0(
      "Incomplete ODE integration: istate=", solver_status[1],
      "; final time=", reached_time, "; requested=", target_time)))
  }
  
  y_hat_i <- as.numeric(
    ode_output[nrow(ode_output), -1]
  )
  
  if (
    length(y_hat_i) != data_list$NY ||
    any(!is.finite(y_hat_i))
  ) {
    return(
      list(
        success = FALSE,
        reason = "Invalid steady-state output"
      )
    )
  }
  
  flux_attempt <- tryCatch(
    {
      list(
        ok = TRUE,
        value = compute_all_fluxes_R(
          y = y_hat_i,
          parms = parms_i,
          C = parms_i$C,
          acetate = parms_i$acetate
        )
      )
    },
    error = function(e) {
      list(
        ok = FALSE,
        message = conditionMessage(e)
      )
    }
  )
  
  if (!isTRUE(flux_attempt$ok)) {
    return(
      list(
        success = FALSE,
        reason = paste0(
          "Flux calculation failure: ",
          flux_attempt$message
        )
      )
    )
  }
  
  f_hat_i <- flux_attempt$value
  
  if (any(!is.finite(f_hat_i))) {
    return(
      list(
        success = FALSE,
        reason = "Invalid flux output: non-finite value"
      )
    )
  }
  
  # Concentration replication retains the fitted Stan convention: the first
  # observation time for y_rep, and the final time for the flux prediction.
  y_for_rep <- as.numeric(ode_output[2L, -1L])
  y_rep_i <- 10^rnorm(length(y_for_rep), log10(pmax(y_for_rep, 1e-12)), sigma_c_i)
  # Match the fitted log10-scale observation model for the displayed fluxes.
  # As in the historical Figure 14 pipeline, floor latent predictions at 1e-12.
  # Only indices 3:13 are observed/plotted; redox output 14 is not assessed here.
  f_rep_i <- 10^rnorm(length(f_hat_i),
                     mean = log10(pmax(f_hat_i, 1e-12)), sd = sigma_f_i)
  if (any(!is.finite(y_rep_i)) || any(!is.finite(f_rep_i))) {
    return(list(success = FALSE, reason = "Non-finite generated predictive draw"))
  }

  list(
    success = TRUE,
    y_hat = y_hat_i,
    f_hat = f_hat_i,
    y_rep = y_rep_i,
    f_rep = f_rep_i,
    sigma_c = sigma_c_i,
    sigma_f = sigma_f_i,
    theta = theta_i
  )
}


# ============================================================
# 5. Run all prior draws
# ============================================================

# Fail once for missing/wrongly named inputs, before attempting all draws.
check_params <- base_params
check_theta <- raw_draw_to_theta(as.list(setNames(prior_spec$mean, prior_spec$variable)))
check_params[theta_parameter_names] <- as.list(check_theta)
stopifnot(length(ode_func(0, y0_prior, check_params)[[1]]) == 13L,
          all(is.finite(ode_func(0, y0_prior, check_params)[[1]])))

prior_results <- lapply(seq_len(n_prior), function(i) {
  if (i == 1L || i %% 100L == 0L || i == n_prior) {
    cat("Prior simulation", i, "of", n_prior, "\n")
  }
  run_one_prior_draw(i)
})

successful_draws <- vapply(
  prior_results,
  function(x) {
    isTRUE(x$success)
  },
  logical(1)
)

cat(
  "Successful prior simulations:",
  sum(successful_draws),
  "of",
  n_prior,
  "\n"
)

cat(
  "Success rate:",
  round(
    100 * mean(successful_draws),
    2
  ),
  "%\n"
)


# Show the actual reason for failed draws.  The earlier version replaced every
# error by NULL, which hid missing parameters and other coding errors.
if (any(!successful_draws)) {
  
  failure_reasons <- vapply(
    prior_results[!successful_draws],
    function(x) x$reason,
    character(1)
  )
  
  cat("\nFailure reasons:\n")
  print(
    sort(
      table(failure_reasons),
      decreasing = TRUE
    )
  )
}

prior_output_dir <- getOption(
  "tca.prior.output_dir",
  file.path(dirname(stan_model_file), "results",
            paste0("prior_predictive_", format(Sys.time(), "%Y%m%d_%H%M%S")))
)
dir.create(prior_output_dir, recursive = TRUE, showWarnings = FALSE)
prior_output_dir <- normalizePath(prior_output_dir, mustWork = TRUE)
saveRDS(list(prior_draws = prior_draws, prior_results = prior_results,
             successful_draws = successful_draws, prior_spec = prior_spec,
             data_list = data_list, stan_model_file = stan_model_file,
             stan_model_md5 = unname(tools::md5sum(stan_model_file))),
        file.path(prior_output_dir, "prior_simulations.rds"))
write.csv(prior_spec, file.path(prior_output_dir, "prior_spec.csv"), row.names = FALSE)
cat("Prior results saved to:", prior_output_dir, "\n")
if (any(!successful_draws)) {
  warning("Predictive plots use only successful integrations; inspect failure reasons and the saved success flags.",
          call. = FALSE)
}

if (!any(successful_draws)) {
  stop(
    paste0(
      "No prior simulation succeeded. ",
      "Read the printed failure reason above; it now contains the original error."
    ),
    call. = FALSE
  )
}



prior_y_rep <- do.call(
  rbind,
  lapply(
    prior_results[successful_draws],
    `[[`,
    "y_rep"
  )
)

prior_f_rep_all <- do.call(
  rbind,
  lapply(
    prior_results[successful_draws],
    `[[`,
    "f_rep"
  )
)

colnames(prior_y_rep) <- paste0(
  "y_rep[",
  seq_len(data_list$NY),
  "]"
)

colnames(prior_f_rep_all) <- paste0(
  "f_rep[",
  1:14,
  "]"
)


# ============================================================
# 6. Existing posterior predictive draws
# ============================================================

posterior_y_rep <- posterior::as_draws_matrix(
  tca_ferm_fit_real$draws(
    variables = "y_rep"
  )
)

# Reconstruct BOTH predictive distributions from latent fluxes. The base
# fit's cached f_rep uses additive Normal noise, which differs from the
# fitted likelihood; the historical run_fast_stan_prior_predictive.R already
# bypassed it in the same way. No ODE simulation or posterior MCMC is needed.
generate_log10_flux_predictions <- function(f_hat, sigma_f, seed) {
  f_hat <- as.matrix(f_hat)
  stopifnot(ncol(f_hat) == 14L, length(sigma_f) == nrow(f_hat),
            all(is.finite(f_hat)), all(is.finite(sigma_f)), all(sigma_f > 0))
  set.seed(seed)
  predictions <- vapply(seq_len(ncol(f_hat)), function(j) {
    10^rnorm(nrow(f_hat), mean = log10(pmax(f_hat[, j], 1e-12)), sd = sigma_f)
  }, numeric(nrow(f_hat)))
  colnames(predictions) <- paste0("f_rep[", seq_len(ncol(f_hat)), "]")
  predictions
}

# Replot helpers supply the original records from prior_simulations.rds.
# Never take f_rep from an older cache without checking its observation model.
if (!exists("prior_results", inherits = FALSE) ||
    !exists("successful_draws", inherits = FALSE)) {
  if (exists("saved", inherits = FALSE) &&
      all(c("prior_results", "successful_draws") %in% names(saved))) {
    prior_results <- saved$prior_results
    successful_draws <- saved$successful_draws
  } else {
    stop("Load prior_simulations.rds to recover prior_results and successful_draws before replotting.",
         call. = FALSE)
  }
}
prior_records <- prior_results[successful_draws]
prior_f_hat <- do.call(rbind, lapply(prior_records, `[[`, "f_hat"))
prior_sigma_f <- vapply(prior_records, `[[`, numeric(1), "sigma_f")
prior_f_rep_all <- generate_log10_flux_predictions(prior_f_hat, prior_sigma_f,
                                                  seed = 20260817L)

post_flux_draws <- posterior::as_draws_matrix(
  tca_ferm_fit_real$draws(variables = c("f_hat", "sigma_f"))
)
posterior_f_hat <- as.matrix(post_flux_draws[
  , paste0("f_hat[", 1:14, "]"), drop = FALSE
])
posterior_sigma_f <- as.numeric(post_flux_draws[, "sigma_f"])
posterior_f_rep_all <- generate_log10_flux_predictions(posterior_f_hat,
                                                       posterior_sigma_f,
                                                       seed = 20260921L)

# A log10 density has a different unit from an original-scale density. Keep
# all positive predictive samples, including values near the likelihood floor;
# do not reproduce the historical plot's silent value > 1e-10 filtering.
flux_prediction_floor_diagnostics <- rbind(
  data.frame(Distribution = "Prior predictive", flux_index = 1:14,
             plotted = 1:14 %in% 3:13, n_used = nrow(prior_f_hat),
             latent_floor_fraction = colMeans(prior_f_hat <= 1e-12),
             predictive_le_1e_minus10_fraction = colMeans(prior_f_rep_all <= 1e-10)),
  data.frame(Distribution = "Posterior predictive", flux_index = 1:14,
             plotted = 1:14 %in% 3:13, n_used = nrow(posterior_f_hat),
             latent_floor_fraction = colMeans(posterior_f_hat <= 1e-12),
             predictive_le_1e_minus10_fraction = colMeans(posterior_f_rep_all <= 1e-10))
)
write.csv(flux_prediction_floor_diagnostics,
          file.path(prior_output_dir, "flux_prediction_floor_diagnostics.csv"),
          row.names = FALSE)


# Put indexed generated quantities in numerical order.  This prevents a
# possible alphabetical order such as [1], [10], [11], ..., [2].
order_indexed_draws <- function(draws, variable) {
  
  expected_names <- paste0(
    variable,
    "[",
    seq_len(ncol(draws)),
    "]"
  )
  
  if (all(expected_names %in% colnames(draws))) {
    draws <- draws[, expected_names, drop = FALSE]
  }
  
  draws
}

posterior_y_rep <- order_indexed_draws(
  posterior_y_rep,
  "y_rep"
)

posterior_f_rep_all <- order_indexed_draws(
  posterior_f_rep_all,
  "f_rep"
)

# ============================================================
# 7. Summarise predictive distributions
# ============================================================

# Extract the latest observation without assuming whether time is stored in
# rows or columns.  This is the assumption that caused the previous f_obs
# length error.
extract_latest_vector <- function(
    x,
    accepted_sizes,
    n_time = NULL,
    object_name = "object"
) {
  
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  
  dims <- dim(x)
  
  if (is.null(dims)) {
    
    values <- as.numeric(x)
    
    if (!(length(values) %in% accepted_sizes)) {
      stop(
        object_name,
        " is a vector of length ",
        length(values),
        "; expected one of: ",
        paste(accepted_sizes, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    
    return(values)
  }
  
  if (length(dims) != 2L) {
    stop(
      object_name,
      " must be a vector or a two-dimensional matrix/data frame; dimensions are ",
      paste(dims, collapse = " x "),
      ".",
      call. = FALSE
    )
  }
  
  nr <- dims[1]
  nc <- dims[2]
  
  # Prefer the orientation that agrees with N_t when it is informative.
  if (
    !is.null(n_time) &&
    nr == n_time &&
    nc %in% accepted_sizes
  ) {
    return(as.numeric(x[nr, , drop = TRUE]))
  }
  
  if (
    !is.null(n_time) &&
    nc == n_time &&
    nr %in% accepted_sizes
  ) {
    return(as.numeric(x[, nc, drop = TRUE]))
  }
  
  # Fallback when there is only one plausible variable dimension.
  if (nc %in% accepted_sizes && !(nr %in% accepted_sizes)) {
    return(as.numeric(x[nr, , drop = TRUE]))
  }
  
  if (nr %in% accepted_sizes && !(nc %in% accepted_sizes)) {
    return(as.numeric(x[, nc, drop = TRUE]))
  }
  
  # A one-row or one-column object is also unambiguous.
  if (nr == 1L && nc %in% accepted_sizes) {
    return(as.numeric(x[1, , drop = TRUE]))
  }
  
  if (nc == 1L && nr %in% accepted_sizes) {
    return(as.numeric(x[, 1, drop = TRUE]))
  }
  
  stop(
    "Cannot determine the orientation of ",
    object_name,
    " with dimensions ",
    paste(dims, collapse = " x "),
    ". Expected the variable dimension to have size ",
    paste(accepted_sizes, collapse = " or "),
    ".",
    call. = FALSE
  )
}


# Observed concentrations at the final observation time
observed_y <- extract_latest_vector(
  data_list$y_obs,
  accepted_sizes = data_list$NY,
  n_time = data_list$N_t,
  object_name = "data_list$y_obs"
)

stopifnot(
  length(observed_y) == data_list$NY,
  length(observed_y) == ncol(prior_y_rep),
  length(observed_y) == ncol(posterior_y_rep)
)



summarise_predictive <- function(
    draws,
    observed,
    labels,
    distribution
) {
  
  relative_draws <- sweep(
    draws,
    MARGIN = 2,
    STATS = observed,
    FUN = "/"
  )
  
  quantiles <- t(
    apply(
      relative_draws,
      MARGIN = 2,
      quantile,
      probs = c(
        0.05,
        0.25,
        0.5,
        0.75,
        0.95
      ),
      na.rm = TRUE
    )
  )
  
  data.frame(
    Variable = labels,
    Distribution = distribution,
    Lower90 = quantiles[, 1],
    Lower50 = quantiles[, 2],
    Median = quantiles[, 3],
    Upper50 = quantiles[, 4],
    Upper90 = quantiles[, 5]
  )
}


if (data_list$NY == 13) {
  
  concentration_labels <- c(
    "Pyruvate",
    "Oxaloacetate",
    "Acetyl-CoA",
    "Citrate",
    "cis-Aconitate",
    "Isocitrate",
    "\u03b1-Ketoglutarate",
    "Succinyl-CoA",
    "Succinate",
    "Fumarate",
    "Acetyl-phosphate",
    "NAD\u207A",
    "NADH"
  )
  
} else {
  
  concentration_labels <- paste0(
    "Concentration ",
    seq_len(data_list$NY)
  )
}

flux_labels <- c(
  "vCS",
  "vACO1",
  "vACO2",
  "vIDH",
  "vAKGDH",
  "vSCS",
  "vSDH",
  "vFUM",
  "vMDH",
  "vPTA",
  "vACKA"
)



concentration_summary <- rbind(
  summarise_predictive(
    prior_y_rep,
    observed_y,
    concentration_labels,
    "Prior predictive"
  ),
  summarise_predictive(
    posterior_y_rep,
    observed_y,
    concentration_labels,
    "Posterior predictive"
  )
)

# Extract the 11 observed model fluxes:
# CS, ACO1, ACO2, IDH, AKGDH, SCS, SDH, FUM, MDH, PTA, ACKA

prior_f_rep <- prior_f_rep_all[
  ,
  3:13,
  drop = FALSE
]

# Posterior f_rep may contain either the 11 observed fluxes, the 13 reaction
# fluxes (PK through ACKA), or all 14 outputs including the redox flux.
if (ncol(posterior_f_rep_all) == 11) {
  
  posterior_f_rep <- posterior_f_rep_all
  
} else if (ncol(posterior_f_rep_all) %in% c(13, 14)) {
  
  posterior_f_rep <- posterior_f_rep_all[
    ,
    3:13,
    drop = FALSE
  ]
  
} else {
  
  stop(
    "Unexpected number of columns in posterior_f_rep_all: ",
    ncol(posterior_f_rep_all)
  )
}

# Extract the latest observed flux vector.  Depending on how data_list was
observed_f_tca <- as.numeric(
  data_list$f_obs_tca
)

observed_f_ferm <- as.numeric(
  data_list$f_obs_ferm
)

stopifnot(
  length(observed_f_tca) == 1L,
  length(observed_f_ferm) == 1L
)

observed_f <- c(
  rep(observed_f_tca, 9),
  rep(observed_f_ferm, 2)
)

names(observed_f) <- flux_labels

stopifnot(
  ncol(prior_f_rep) == 11L,
  ncol(posterior_f_rep) == 11L,
  length(observed_f) == 11L
)



flux_summary <- rbind(
  summarise_predictive(
    prior_f_rep,
    observed_f,
    flux_labels,
    "Prior predictive"
  ),
  summarise_predictive(
    posterior_f_rep,
    observed_f,
    flux_labels,
    "Posterior predictive"
  )
)

concentration_summary$Variable <- factor(
  concentration_summary$Variable,
  levels = rev(concentration_labels)
)

flux_summary$Variable <- factor(
  flux_summary$Variable,
  levels = rev(flux_labels)
)



library(ggplot2)
library(dplyr)

prior_fill <- "grey65"
posterior_fill <- "#56B4E9"
ci_colour <- "#005A9C"


# ============================================================
# Convert a predictive matrix to long format
# ============================================================

predictive_matrix_to_long <- function(
    draws,
    labels,
    distribution
) {
  
  # Drop the posterior::draws_matrix class completely.  Otherwise the
  # value column can retain that class and fail when rows are combined.
  draws <- matrix(
    as.numeric(draws),
    nrow = nrow(draws),
    ncol = ncol(draws),
    dimnames = dimnames(draws)
  )
  
  stopifnot(
    ncol(draws) == length(labels)
  )
  
  data.frame(
    Variable = factor(
      rep(
        labels,
        each = nrow(draws)
      ),
      levels = labels
    ),
    Distribution = distribution,
    value = as.numeric(draws)
  )
}


# ============================================================
# Predictive-density plotting function
# ============================================================

make_predictive_density_plot <- function(
    prior_draws,
    posterior_draws,
    observed,
    labels,
    x_label,
    ncol = 5,
    x_limits = NULL,
    facet_scales = "free_x",
    interval_prob = 0.90,
    display_prob = 0.95,
    log_scale = TRUE
) {
  
  prior_draws <- matrix(
    as.numeric(prior_draws),
    nrow = nrow(prior_draws),
    ncol = ncol(prior_draws),
    dimnames = dimnames(prior_draws)
  )
  
  posterior_draws <- matrix(
    as.numeric(posterior_draws),
    nrow = nrow(posterior_draws),
    ncol = ncol(posterior_draws),
    dimnames = dimnames(posterior_draws)
  )
  observed <- as.numeric(observed)
  
  stopifnot(
    ncol(prior_draws) == length(labels),
    ncol(posterior_draws) == length(labels),
    length(observed) == length(labels),
    all(observed > 0)
  )
  
  
  # ----------------------------------------------------------
  # Long-format predictive draws
  # ----------------------------------------------------------
  
  df_all <- rbind(
    predictive_matrix_to_long(
      prior_draws,
      labels,
      "Prior predictive"
    ),
    predictive_matrix_to_long(
      posterior_draws,
      labels,
      "Posterior predictive"
    )
  )
  
  # Both concentration and flux figures use log10 predictions and coordinates.
  df_all <- df_all %>% filter(is.finite(value))
  if (log_scale) df_all <- df_all %>% filter(value > 0)
  df_all <- df_all %>% mutate(plot_value = if (log_scale) log10(value) else value)

  # Quantiles set a view window only. Every valid draw remains available
  # for the KDE and posterior interval; displayed densities are not renormalized.
  stopifnot(display_prob > 0, display_prob <= 1,
            interval_prob > 0, interval_prob < 1)
  display_tail <- (1 - display_prob) / 2
  display_bounds <- df_all %>%
    group_by(Variable, Distribution) %>%
    summarise(
      display_lower = quantile(plot_value, display_tail, na.rm = TRUE),
      display_upper = quantile(plot_value, 1 - display_tail, na.rm = TRUE),
      .groups = "drop"
    )

  # ----------------------------------------------------------
  # Observed values
  # ----------------------------------------------------------
  
  df_observed <- data.frame(
    Variable = factor(
      labels,
      levels = labels
    ),
    observed_value = observed,
    observed_plot = if (log_scale) log10(observed) else observed
  )
  
  
  # ----------------------------------------------------------
  # Posterior predictive median and central PI
  # ----------------------------------------------------------
  
  posterior_for_pi <- posterior_draws
  posterior_for_pi[
    !is.finite(posterior_for_pi) |
      (log_scale & posterior_for_pi <= 0)
  ] <- NA_real_
  
  posterior_plot <- if (log_scale) log10(posterior_for_pi) else posterior_for_pi
  
  interval_tail <- (1 - interval_prob) / 2
  
  posterior_quantiles <- t(
    apply(
      posterior_plot,
      MARGIN = 2,
      quantile,
      probs = c(
        interval_tail,
        0.5,
        1 - interval_tail
      ),
      na.rm = TRUE
    )
  )
  
  df_pi <- data.frame(
    Variable = factor(
      labels,
      levels = labels
    ),
    pi_lower = posterior_quantiles[, 1],
    median_value = posterior_quantiles[, 2],
    pi_upper = posterior_quantiles[, 3]
  )
  
  
  # Each facet includes the central range of BOTH distributions, the
  # observed value, and the entire posterior PI. The latter remains 90%.
  view_ranges <- display_bounds %>%
    group_by(Variable) %>%
    summarise(view_lower = min(display_lower),
              view_upper = max(display_upper), .groups = "drop") %>%
    left_join(df_observed[, c("Variable", "observed_plot")], by = "Variable") %>%
    left_join(df_pi[, c("Variable", "pi_lower", "pi_upper")], by = "Variable") %>%
    mutate(view_lower = pmin(view_lower, observed_plot, pi_lower),
           view_upper = pmax(view_upper, observed_plot, pi_upper))
  if (!is.null(x_limits)) {
    stopifnot(length(x_limits) == 2L, all(is.finite(x_limits)), diff(x_limits) > 0)
    view_ranges$view_lower <- x_limits[1]
    view_ranges$view_upper <- x_limits[2]
  }

  # Precompute KDEs from all draws and evaluate them over each panel's view.
  # Supplying a shorter evaluation grid does not discard samples or change
  # their weights/bandwidth. This supports independent facets without extras.
  density_groups <- split(df_all, interaction(df_all$Variable,
                                               df_all$Distribution, drop = TRUE))
  density_metadata <- vector("list", length(density_groups))
  density_frames <- vector("list", length(density_groups))
  for (i in seq_along(density_groups)) {
    group <- density_groups[[i]]
    values <- group$plot_value
    bounds <- view_ranges[match(as.character(group$Variable[1]),
                                as.character(view_ranges$Variable)), ]
    if (length(values) < 2L) {
      stop("At least two valid predictive draws are required per density.", call. = FALSE)
    }
    kde <- stats::density(values, bw = "nrd0", adjust = 1, n = 1024L,
                          from = bounds$view_lower, to = bounds$view_upper)
    density_frames[[i]] <- data.frame(
      Variable = factor(as.character(group$Variable[1]), levels = labels),
      Distribution = group$Distribution[1],
      plot_value = kde$x, density = kde$y
    )
    density_metadata[[i]] <- data.frame(
      Variable = as.character(group$Variable[1]),
      Distribution = group$Distribution[1],
      n_used = length(values), bandwidth = kde$bw,
      view_lower = bounds$view_lower, view_upper = bounds$view_upper,
      outside_fraction = mean(values < bounds$view_lower | values > bounds$view_upper)
    )
  }
  df_density <- do.call(rbind, density_frames)
  density_metadata <- do.call(rbind, density_metadata)

  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  pi_legend_label <- paste0(
    "Posterior ",
    round(100 * interval_prob),
    "% PI"
  )
  
  plot <- ggplot() +
    
    # Light-blue posterior predictive interval band
    geom_rect(
      data = df_pi,
      mapping = aes(
        xmin = pi_lower,
        xmax = pi_upper,
        ymin = -Inf,
        ymax = Inf,
        fill = pi_legend_label
      ),
      alpha = 0.32,
      colour = NA,
      inherit.aes = FALSE
    ) +
    
    # Prior predictive density
    geom_area(
      data = df_density %>%
        filter(
          Distribution == "Prior predictive"
        ),
      mapping = aes(
        x = plot_value,
        y = density,
        fill = Distribution
      ),
      colour = prior_fill,
      linewidth = 0.35,
      alpha = 0.30,
      stat = "identity",
      position = "identity",
      na.rm = TRUE
    ) +
    
    # Posterior predictive density
    geom_area(
      data = df_density %>%
        filter(
          Distribution == "Posterior predictive"
        ),
      mapping = aes(
        x = plot_value,
        y = density,
        fill = Distribution
      ),
      colour = posterior_fill,
      linewidth = 0.45,
      alpha = 0.45,
      stat = "identity",
      position = "identity",
      na.rm = TRUE
    ) +
    
    # Posterior predictive median: blue dashed vertical line
    geom_vline(
      data = df_pi,
      mapping = aes(
        xintercept = median_value,
        colour = "Posterior median",
        linetype = "Posterior median"
      ),
      linewidth = 0.65
    ) +
    
    # Observed value: orange solid vertical line
    geom_vline(
      data = df_observed,
      mapping = aes(
        xintercept = observed_plot,
        colour = "Observed",
        linetype = "Observed"
      ),
      linewidth = 0.65
    ) +
    
    facet_wrap(
      facets = vars(Variable),
      ncol = ncol,
      scales = facet_scales
    ) +
    
    scale_fill_manual(
      name = NULL,
      values = c(
        "Prior predictive" = prior_fill,
        "Posterior predictive" = posterior_fill,
        "Posterior 90% PI" = "#D9EEF8"
      )
    ) +
    
    scale_colour_manual(
      name = NULL,
      values = c(
        "Posterior median" = ci_colour,
        "Observed" = "#E66101"
      )
    ) +
    
    scale_linetype_manual(
      name = NULL,
      values = c(
        "Posterior median" = "dashed",
        "Observed" = "solid"
      )
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.07, 0.07))
    ) +
    
    scale_y_continuous(
      breaks = scales::breaks_pretty(
        n = 3
      ),
      expand = expansion(
        mult = c(0.03, 0.08)
      )
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      clip = "on"
    ) +
    
    labs(
      x = x_label,
      y = "Predictive density",
      caption = paste0("Central ", round(100 * display_prob),
                       "% view; tails outside each panel are not shown. Posterior PI: ",
                       round(100 * interval_prob), "%.")
    ) +
    
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(
          alpha = c(0.30, 0.45, 0.32)
        )
      ),
      colour = guide_legend(
        order = 2,
        override.aes = list(linewidth = 0.8)
      ),
      linetype = guide_legend(
        order = 2
      )
    ) +
    
    theme_classic(
      base_size = 9,
      base_family = "sans"
    ) +
    
    theme(
      plot.caption = element_text(size = 7, hjust = 0),
      legend.position = "top",
      legend.justification = "center",
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        size = 9,
        face = "plain",
        colour = "black"
      ),
      
      axis.title = element_text(
        size = 9,
        colour = "black"
      ),
      
      axis.text = element_text(
        size = 7,
        colour = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.35,
        colour = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.4,
        colour = "black"
      ),
      
      panel.spacing = grid::unit(
        1,
        "lines"
      ),
      
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
  attr(plot, "predictive_view") <- list(
    display_prob = display_prob, interval_prob = interval_prob,
    ranges = view_ranges, density_metadata = density_metadata,
    posterior_intervals = df_pi
  )
  plot

}



concentration_density_plot <- make_predictive_density_plot(
  prior_draws = prior_y_rep,
  posterior_draws = posterior_y_rep,
  observed = observed_y,
  labels = concentration_labels,
  x_label = expression(
    log[10] * " predicted concentration (mM)"
  ),
  ncol = 5,
  x_limits = NULL,
  facet_scales = "free_x",
  interval_prob = 0.90,
  display_prob = 0.95
)



flux_density_plot <- make_predictive_density_plot(
  prior_draws = prior_f_rep,
  posterior_draws = posterior_f_rep,
  observed = observed_f,
  labels = flux_labels,
  x_label = expression(
    log[10] * " flux (mM " * s^{-1} * ")"
  ),
  ncol = 5,
  x_limits = NULL,
  facet_scales = "free_x",
  interval_prob = 0.90,
  display_prob = 0.97,
  log_scale = TRUE
)

print(concentration_density_plot)
print(flux_density_plot)


library(patchwork)

final_fig <- (
  (concentration_density_plot + labs(caption = NULL)) /
    (flux_density_plot + labs(caption = NULL))
) +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  )

final_fig <- final_fig &
  theme(
    legend.position = "top",
    
    plot.tag = element_text(
      size = 10,
      face = "bold",
      colour = "black"
    ),
    
    plot.tag.position = c(0, 1)
  )

print(final_fig)


# Save data/plot objects first, so a later export never loses the simulation.
saveRDS(list(prior_y_rep = prior_y_rep, prior_f_rep_all = prior_f_rep_all,
             posterior_y_rep = posterior_y_rep, posterior_f_rep_all = posterior_f_rep_all,
             flux_observation_model = "Normal on log10 flux; latent floor 1e-12",
             prior_f_hat = prior_f_hat, prior_sigma_f = prior_sigma_f,
             posterior_f_hat = posterior_f_hat, posterior_sigma_f = posterior_sigma_f,
             flux_prediction_floor_diagnostics = flux_prediction_floor_diagnostics,
             concentration_summary = concentration_summary, flux_summary = flux_summary),
        file.path(prior_output_dir, "predictive_draws_and_summaries.rds"))
saveRDS(list(concentration = concentration_density_plot, flux = flux_density_plot,
             combined = final_fig), file.path(prior_output_dir, "predictive_plots.rds"))

save_prior_plot <- function(plot, filename, width, height) {
  ggplot2::ggsave(file.path(prior_output_dir, paste0(filename, ".png")),
                 plot = plot, width = width, height = height, units = "cm",
                 dpi = 300, bg = "white")
}
save_prior_plot(concentration_density_plot, "priorposterior_concentrations1024", 20, 10)
save_prior_plot(flux_density_plot, "priorposterior_fluxes1024", 20, 10)
save_prior_plot(final_fig, "priorposterior_combined1024", 20, 20)
cat("Predictive draws and figures saved to:", prior_output_dir, "\n")

