library(cmdstanr)
install_cmdstan()   
library(jsonlite)

model_real <- cmdstan_model("tcaferm_real.stan")


y_obs <- c(
  pyruvate = 3.66,
  oxaloacetate = 0.000487,
  acetyl_coa = 0.606,
  citrate = 1.96,
  cis_aconitate = 0.0161,
  isocitrate = 0.0367,
  aketo = 0.443,
  succinyl_coa = 0.233,
  succinate = 0.569,
  fumarate = 0.115,
  acetyl_p = 1.07,
  nad = 2.8,
  nadh = 0.0545
)

f_obs_ferm <- 1.0833
f_obs_tca <- 0.3



params <- list(
  ki_nadh_pdh = 0.001,
  ki_nadh_cs = 0.0028,
  
  ki_atp_pk = 4,
  
  
  # michaelis constant
  # pk
  km_pyruvate_pk = 0.0505,
  km_pep_pk = 0.2666,
  km_atp_pk = 16.1752,
  km_adp_pk = 0.1919,
  
  
  # pdh
  km_pyruvate_pdh = 0.3296, 
  km_coa_pdh = 0.0471, 
  km_nad_pdh = 0.0619, 
  km_nadh_pdh = 0.215,
  km_acetyl_coa_pdh = 0.215,
  km_co2_pdh = 0.215,
  
  # pta
  km_acetyl_coa_pta = 0.0491,
  km_pi_pta = 1.9606,
  km_acetyl_p_pta = 0.6297,
  km_coa_pta = 0.0521,
  
  # ack
  km_acetyl_p_ack = 0.0938,
  km_adp_ack = 0.2799,
  km_acetate_ack = 7.9678,
  km_atp_ack = 0.096,
  
  
  # cs
  km_citrate_cs = 0.0018, 
  km_coa_sh_cs = 44914.6131, # 0.0027, 
  km_oxaloacetate_cs = 0.0119, 
  km_acetyl_coa_cs = 0.0521, 
  
  
  # aconitase
  km_citrate_aco = 0.6464, 
  km_cis_aconitate_aco1 = 0.0095,
  km_cis_aconitate_aco2 = 0.2098,
  km_isocitrate_aco = 0.0872, 
  
  # isocitrate dehydrogenase
  km_isocitrate_id = 0.0109, 
  km_nad_id = 1.5211, 
  km_aketo_id = 0.4744, 
  km_nadh_id = 0.007, # 24.3, 
  km_co2_id = 2.4317,
  
  
  
  # akd checked not much different
  km_nad_akd = 0.0512, 
  km_aketo_akd = 0.0026, 
  km_coa_sh_akd = 0.0512, 
  km_nadh_akd = 0.1953, 
  km_succinyl_coa_akd = 0.1953, 
  km_co2_akd = 0.1953,
  
  
  # scs finished check 
  km_coa_sh_scs = 2.0659, 
  km_succinyl_coa_scs = 0.0316, 
  km_adp_scs = 0.3244, 
  km_succinate_scs = 0.1444, 
  km_atp_scs = 0.0426, 
  km_pi_scs = 0.8953,
  
  
  # sd finished check 
  km_succinate_sd = 0.0021,  
  km_fad_sd = 0.002, # q
  km_fumarate_sd = 0.0233, 
  km_fadh2_sd = 12.5427, # qh2
  
  
  # fumarase
  km_fumarate_fumarase = 0.4222, 
  km_malate_fumarase = 0.6905, 
  
  # malate dehydrogenase dont change same ecoli exp
  km_nad_md = 0.418, 
  km_oxaloacetate_md = 0.0244, 
  km_malate_md = 3.8086, 
  km_nadh_md = 0.0301, 
  
  
  
  
  # forward and reverse k_cat 
  
  kf_pk = 555.4713,
  kr_pk = 0.4124,
  kf_pdh = 36.1192, 
  kr_pdh = 0,
  kf_cs = 141.5176, 
  kr_cs = 5.8632, 
  kf1_aco = 5.7836, 
  kf2_aco = 4.5447,
  kr1_aco = 3.8424, 
  kr2_aco = 1.0233,
  kf_id = 82.0188, 
  kr_id = 14.7374,
  kf_akd = 2.0345, 
  kr_akd = 0.0011,
  kf_scs = 35.8092, 
  kr_scs = 60.5031, 
  kf_sd = 32.8092, 
  kr_sd = 180.5162, 
  kf_fumarase = 723.0053, 
  kr_fumarase = 395.3746,
  kf_md = 132.9511, 
  kr_md = 1557.9368, 
  kf_pta = 22.1455, 
  kr_pta = 237.7203,
  kf_ack = 4791.0944,
  kr_ack = 852.687,
  
  
  
  conserved_coa = 2.2090,
  conserved_quantity2 = 5.0533,
  
  
  
  nadh_constant = 0.0545,
  atp_constant = 8.13,
  fadh2_constant = 0.0173,
  co2_constant = 0.0752,
  
  adp_constant = 0.704,
  fad_constant = 0.173,
  nad_constant = 2.80,
  pep_constant = 0.231,
  pi_constant = 23.9,
  
  
  pyruvate_initial = 3.66,
  oxaloacetate_initial = 0.000487,
  acetyl_coa_initial = 0.606,
  citrate_initial = 1.96,
  cis_aconitate_initial = 0.0161,
  isocitrate_initial = 0.0367,
  aketo_initial = 0.443,
  coa_initial = 1.37,
  succinyl_coa_initial = 0.233,
  succinate_initial = 0.569,
  fumarate_initial = 0.115,
  malate_initial = 1.68,
  
  acetyl_p_initial = 1.07
  
  
)



build_stan_data_from_params <- function(params, t_obs) {
  # ==== 1. 定义物种数 ====
  NY <- 13   # 你 Stan 里 RHS 里 dy[1..10]
  
  # ==== 2. 初始条件 y0 ====
  y0 <- c(
    params$pyruvate_initial,
    params$oxaloacetate_initial,
    params$acetyl_coa_initial,
    params$citrate_initial,
    params$cis_aconitate_initial,
    params$isocitrate_initial,
    params$aketo_initial,
    params$succinyl_coa_initial,
    params$succinate_initial,
    params$fumarate_initial,
    params$acetyl_p_initial,
    params$nad_constant,
    params$nadh_constant
  )
  
  # ==== 3. 常量池 ====
  # conserved quantities + cofactors
  C <- c(
    params$nadh_constant,
    params$nad_constant,
    params$atp_constant,
    params$adp_constant,
    params$fad_constant,
    params$fadh2_constant,
    params$co2_constant,
    params$pep_constant,
    params$conserved_coa,
    params$conserved_quantity2,
    params$pi_constant
  )
  
  # ==== 4. 分模块参数打包 ====
  # 按照 Stan 文件里 slice_from 的顺序来排
  p_pdh <- c(
    params$km_co2_pdh, 
    params$km_pyruvate_pdh, 
    params$km_coa_pdh, 
    params$km_nad_pdh, 
    params$km_acetyl_coa_pdh, 
    params$km_nadh_pdh, 
    params$kf_pdh, 
    params$ki_nadh_pdh,
    params$kr_pdh
  )
  
  p_cs <- c(
    params$km_oxaloacetate_cs, 
    params$km_acetyl_coa_cs, 
    params$km_citrate_cs, 
    params$km_coa_sh_cs, 
    params$kf_cs, 
    params$kr_cs, 
    params$ki_nadh_cs
  )
  
  p_aco <- c(
    params$km_citrate_aco, 
    params$km_cis_aconitate_aco1, 
    params$km_cis_aconitate_aco2, 
    params$km_isocitrate_aco,
    params$kf1_aco, 
    params$kr1_aco, 
    params$kf2_aco, 
    params$kr2_aco
  )
  
  p_id <- c(
    params$km_co2_id, 
    params$km_isocitrate_id, 
    params$km_nadh_id, 
    params$km_nad_id, 
    params$km_aketo_id,
    params$kf_id,
    params$kr_id
  )
  
  p_akd <- c(
    params$km_co2_akd, 
    params$km_aketo_akd, 
    params$km_nad_akd, 
    params$km_nadh_akd, 
    params$km_coa_sh_akd, 
    params$km_succinyl_coa_akd,
    params$kf_akd,
    params$kr_akd
  )
  
  p_scs <- c(
    params$km_succinyl_coa_scs, 
    params$km_adp_scs, 
    params$km_succinate_scs, 
    params$km_coa_sh_scs, 
    params$km_atp_scs, 
    params$kf_scs, 
    params$kr_scs,
    params$km_pi_scs
  )
  
  p_sd <- c(
    params$km_succinate_sd, 
    params$km_fad_sd, 
    params$km_fumarate_sd, 
    params$km_fadh2_sd, 
    params$kf_sd, 
    params$kr_sd
  )
  
  p_fum <- c(
    params$km_fumarate_fumarase,
    params$km_malate_fumarase,
    params$kf_fumarase,
    params$kr_fumarase
  )
  
  p_md <- c(
    params$km_nad_md, 
    params$km_oxaloacetate_md,
    params$km_malate_md,
    params$km_nadh_md,
    params$kf_md,
    params$kr_md
  )
  
  p_pk <- c(
    params$km_pep_pk,
    params$km_atp_pk,
    params$km_adp_pk,
    params$kf_pk,
    params$km_pyruvate_pk,
    params$ki_atp_pk,
    params$kr_pk
  )
  
  p_pta <- c(
    params$km_acetyl_coa_pta, 
    params$km_pi_pta, 
    params$km_acetyl_p_pta, 
    params$km_coa_pta, 
    params$kf_pta, 
    params$kr_pta
  )
  
  p_ack <- c(
    params$km_acetyl_p_ack, 
    params$km_adp_ack, 
    params$km_acetate_ack, 
    params$km_atp_ack, 
    params$kf_ack, 
    params$kr_ack
  )
  
  
  
  list(
    NY = NY,
    N_t = length(t_obs),
    t_obs = t_obs,
    y_obs = y_obs,
    f_obs_ferm = f_obs_ferm,
    f_obs_tca = f_obs_tca,
    y0 = y0,
    C = C,
    p_pdh = p_pdh,
    p_cs = p_cs,
    p_aco = p_aco,
    p_id = p_id,
    p_akd = p_akd,
    p_scs = p_scs,
    p_sd = p_sd,
    p_fum = p_fum,
    p_md = p_md,
    p_pk = p_pk,
    p_pta = p_pta,
    p_ack = p_ack
  )
}


T_end <- 60*60*2400

t_obs <- c(as.double(T_end))

data_list <- build_stan_data_from_params(params, t_obs)
data_list$y0 <- as.vector(data_list$y0)  # num [1:10]

# y_obs 必须是二维矩阵 (N_t × NY)
data_list$y_obs <- matrix(data_list$y_obs, nrow = data_list$N_t, ncol = data_list$NY)

# 检查一下
str(data_list$y0)   # num [1:10]
str(data_list$y_obs) # num [1, 1:10]



tca_ferm_fit_real <- model_real$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 8000,
  adapt_delta = 0.9,
  max_treedepth = 12,
  seed = 1234
) 


tca_ferm_fit_real_test <- model_real$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.90,
  max_treedepth = 12,
  seed = 1234
)


print(tca_ferm_fit_real$summary(), n = Inf)



# 给这次 subunit 修改后的结果单独建文件夹
out_dir <- file.path(
  "results",
  paste0("stan_subunit_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(out_dir, recursive = TRUE)

# 1. 保存各条 chain 的原始 CSV，之后做 generated quantities 也能用
tca_ferm_fit_real$save_output_files(dir = out_dir)

# 2. 保存完整 fit：后验采样 + 采样诊断
tca_ferm_fit_real$save_object(
  file = file.path(out_dir, "stan_fit.rds")
)

# 3. 单独保存后验 draws，保留 chain / iteration 信息
saveRDS(
  tca_ferm_fit_real$draws(),
  file.path(out_dir, "stan_draws.rds")
)

# 4. 保存 summary 表
write.csv(
  tca_ferm_fit_real$summary(),
  file.path(out_dir, "stan_summary.csv"),
  row.names = FALSE
)

# 看实际存到了哪里
normalizePath(out_dir)


library(dplyr)

draws_df <- tca_ferm_fit_real$draws(format = "df")

# 如果你只想保留参数，不要 lp__, chain, iteration 这些
theta_draws <- draws_df %>%
  select(-.chain, -.iteration, -.draw)

write.csv(theta_draws, "posterior1024.csv", row.names = FALSE)

