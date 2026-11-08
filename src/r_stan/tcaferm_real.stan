functions {

  real atp_inhibition(
      real atp_con, 
      real ki_atp) {
        return 1.0; // / (1.0 + atp_con / ki_atp);
  }
  
  real nadh_inhibition(
      real nadh_con, 
      real ki_nadh) {
        return 1.0 / (1.0 + nadh_con / ki_nadh);
  }
  
  vector compute_v_aco(
      real citrate_con, real cis_aconitate_con, real isocitrate_con,
      real km_cit, real km_cis1, real km_cis2, real km_iso,
      real enzyme_aco, real kf1, real kr1, real kf2, real kr2) {
  
        real cit_tilde = citrate_con / km_cit;
        real cis_tilde1 = cis_aconitate_con / km_cis1;
        real cis_tilde2 = cis_aconitate_con / km_cis2;
        real iso_tilde = isocitrate_con   / km_iso;
    
        real denom1 = 1.0 + cit_tilde + cis_tilde1;
        real denom2 = 1.0 + cis_tilde2 + iso_tilde;
    
        vector[2] v;
        v[1] = (kf1 * cit_tilde - kr1 * cis_tilde1) / denom1 * enzyme_aco;
        v[2] = (kf2 * cis_tilde2 - kr2 * iso_tilde) / denom2 * enzyme_aco;
        return v;
   }
  
  real compute_v_akd(
      real a_keto_con,          
      real nad_con,             
      real coa_sh_con,          
      real succinyl_coa_con,    
      real nadh_con,            
      real co2_con,             
      real km_co2, real km_ak, real km_nad, real km_nadh,
      real km_coa, real km_succinyl_coa,
      real enzyme_akd,          
      real kf_akd,
      real kr_akd
    ) {
  
        real ak_tilde          = a_keto_con        / km_ak;
        real nad_tilde         = nad_con           / km_nad;
        real coa_tilde         = coa_sh_con        / km_coa;
        real succinyl_coa_tilde= succinyl_coa_con  / km_succinyl_coa;
        real nadh_tilde        = nadh_con          / km_nadh;
        real co2_tilde         = co2_con           / km_co2;
    
    
        real denom = (1.0 + nad_tilde) * (1.0 + ak_tilde) * (1.0 + coa_tilde)
                   + (1.0 + succinyl_coa_tilde) * (1.0 + nadh_tilde) * (1.0 + co2_tilde)
                   - 1.0;
                   
        real fwd = kf_akd * ak_tilde * coa_tilde * nad_tilde;
        real rev = kr_akd * succinyl_coa_tilde * nadh_tilde * co2_tilde;
    
        real v = (fwd - rev) / denom * enzyme_akd;
    
        return v;
  }
  
  real compute_v_cs(
      real oxaloacetate_con, real acetyl_coa_con, real citrate_con, real coa_sh_con,
      real nadh_con,
      real km_ox, real km_acc, real km_cit, real km_coa_sh,
      real enzyme_cs, real kf_cs, real kr_cs,
      real ki_nadh) {
  
        real ox_tilde  = oxaloacetate_con / km_ox;
        real acc_tilde = acetyl_coa_con   / km_acc;     
        real cit_tilde = citrate_con      / km_cit;
        real coa_tilde = coa_sh_con       / km_coa_sh;
    
        real denom = (1.0 + ox_tilde) * (1.0 + acc_tilde)
                   + (1.0 + cit_tilde) * (1.0 + coa_tilde) - 1.0;
    
    
        real fwd = kf_cs * ox_tilde * acc_tilde;
        real rev = kr_cs * cit_tilde * coa_tilde;
    
    
        real v = (fwd - rev) / denom * enzyme_cs;
        v *= nadh_inhibition(nadh_con, ki_nadh);
  
        return v;
   }
  
  real compute_v_fumarase(
      real fumarate_con, real malate_con,
      real km_fum, real km_mal,
      real enzyme_fum, real kf_fum, real kr_fum) {
  
        real fum_tilde = fumarate_con / km_fum;
        real mal_tilde = malate_con   / km_mal;
    
        real denom = 1.0 + fum_tilde + mal_tilde;   
    
        real v = (kf_fum * fum_tilde - kr_fum * mal_tilde) / denom * enzyme_fum;
  
        return v;
  }
  
  real compute_v_id(
      real isocitrate_con, real nad_con, real nadh_con,
      real aketo_con, real co2_con,
      real km_co2, real km_iso, real km_nadh, real km_nad, real km_aketo,
      real enzyme_id, real kf_id, real kr_id) {
  
        real iso_tilde   = isocitrate_con / km_iso;
        real nad_tilde   = nad_con        / km_nad;
        real nadh_tilde  = nadh_con       / km_nadh;
        real aketo_tilde = aketo_con      / km_aketo;
        real co2_tilde   = co2_con        / km_co2;
    
     
        real denom = (1.0 + iso_tilde) * (1.0 + nad_tilde)
                   + (1.0 + nadh_tilde) * (1.0 + aketo_tilde) * (1.0 + co2_tilde)
                   - 1.0;
    
        real fwd = kf_id * iso_tilde * nad_tilde;
        real rev = kr_id * nadh_tilde * aketo_tilde * co2_tilde;
        
        
        
        real v = (fwd - rev) / denom * enzyme_id;
    
        return v;
  }
  
  real compute_v_md(
      real malate_con, real nad_con, real oaa_con, real nadh_con,
      real km_mal, real km_nad, real km_oaa, real km_nadh,
      real enzyme_md, real kf_md, real kr_md) {
  
        real mal_tilde  = malate_con / km_mal;
        real nad_tilde  = nad_con    / km_nad;
        real oaa_tilde  = oaa_con    / km_oaa;
        real nadh_tilde = nadh_con   / km_nadh;
  
        real denom = (1.0 + mal_tilde) * (1.0 + nad_tilde)
                   + (1.0 + oaa_tilde) * (1.0 + nadh_tilde) - 1.0;
       
        real v = (kf_md * mal_tilde * nad_tilde - kr_md * oaa_tilde * nadh_tilde) / denom * enzyme_md;
    
        return v;
  }
  
  real compute_v_pdh(
      real pyruvate_con, real coa_con, real nad_con,
      real acetyl_coa_con, real nadh_con, real co2_con,
      real km_co2, real km_pyruvate, real km_coa, real km_nad,
      real km_acetyl_coa, real km_nadh,
      real enzyme_pdh, real kf_pdh, real kr_pdh, real ki_nadh) {
  
        real pyruvate_tilde   = pyruvate_con    / km_pyruvate;
        real coa_tilde        = coa_con         / km_coa;
        real nad_tilde        = nad_con         / km_nad;
        real acetyl_coa_tilde = acetyl_coa_con  / km_acetyl_coa;
        real nadh_tilde       = nadh_con        / km_nadh;
        real co2_tilde        = co2_con         / km_co2;
    
        real denom = (1.0 + pyruvate_tilde) * (1.0 + coa_tilde) * (1.0 + nad_tilde)
                   + (1.0 + acetyl_coa_tilde) * (1.0 + nadh_tilde) * (1.0 + co2_tilde)
                   - 1.0;
        
        real fwd = kf_pdh * pyruvate_tilde * coa_tilde * nad_tilde;
        real rev = kr_pdh * acetyl_coa_tilde * nadh_tilde * co2_tilde;
  
        real v = (fwd - rev) / denom * enzyme_pdh;
        v *= nadh_inhibition(nadh_con, ki_nadh);
    
        return v;
  }
  
  real compute_v_pk(
      real km_pyruvate, real km_atp, real km_pep, real km_adp,
      real kf_pk, real kr_pk, real enzyme_pk,
      real pyruvate_con, real atp_con, real pep_con, real adp_con, real ki_atp) {
  
        real pyruvate_tilde = pyruvate_con / km_pyruvate;
        real pep_tilde      = pep_con      / km_pep;
        real atp_tilde      = atp_con      / km_atp;
        real adp_tilde      = adp_con      / km_adp;
    
        real denom = (1.0 + pyruvate_tilde) * (1.0 + atp_tilde)
                   + (1.0 + pep_tilde)      * (1.0 + adp_tilde) - 1.0;
    
        real v = (kf_pk * pep_tilde * adp_tilde - kr_pk * pyruvate_tilde * atp_tilde) / denom * enzyme_pk;
        
        // v*= atp_inhibition(atp_con, ki_atp);
        
        return v;
  }
  
  
  real compute_v_scs(
      real succinyl_coa_con, real adp_con, real succinate_con,
      real coa_sh_con, real nadh_con, real atp_con, real pi_con,
      real km_succinyl, real km_adp, real km_succinate,
      real km_coa, real km_atp, real km_pi,
      real enzyme_scs, real kf_scs, real kr_scs) {
  
        real sucCoA_tilde = succinyl_coa_con / km_succinyl;
        real adp_tilde    = adp_con          / km_adp;
        real suc_tilde    = succinate_con    / km_succinate;
        real coa_tilde    = coa_sh_con       / km_coa;
        real atp_tilde    = atp_con          / km_atp;
        real pi_tilde     = pi_con           / km_pi;
        
    
        real denom = (1.0 + sucCoA_tilde) * (1.0 + adp_tilde) * (1.0 + pi_tilde)
                   + (1.0 + suc_tilde)   * (1.0 + coa_tilde) * (1.0 + atp_tilde)
                   - 1.0;
        
        real v = (kf_scs * sucCoA_tilde * adp_tilde * pi_tilde
                 - kr_scs * suc_tilde * coa_tilde * atp_tilde) / denom * enzyme_scs;
    
        return v;
  }
  
  real compute_v_sd(
      real succinate_con, real uq_con, real fumarate_con, real uqh2_con,
      real km_suc, real km_uq, real km_fum, real km_uqh2,
      real enzyme_sd, real kf_sd, real kr_sd) {
  
        real suc_tilde  = succinate_con / km_suc;
        real uq_tilde   = uq_con       / km_uq;
        real fum_tilde  = fumarate_con / km_fum;
        real uqh2_tilde = uqh2_con     / km_uqh2;
    
        real denom = (1.0 + suc_tilde) * (1.0 + uq_tilde) + (1.0 + fum_tilde) * (1.0 + uqh2_tilde) - 1.0;
        
        real v = (kf_sd * suc_tilde * uq_tilde
                 - kr_sd * fum_tilde * uqh2_tilde) / denom * enzyme_sd;
    
        return v;
  }
  
  real compute_v_pta(real acetyl_coa_con, real pi_con, real acetyl_p_con, real coa_sh_con, 
                          real km_acetyl_coa_pta, real km_pi_pta, real km_acetyl_p_pta, real km_coa_pta, 
                          real kf_pta, real kr_pta, real enzyme_pta) {
                            
    real acetyl_coa_tilde = acetyl_coa_con / km_acetyl_coa_pta;
    real pi_tilde = pi_con                 / km_pi_pta;
    real acetyl_p_tilde = acetyl_p_con     / km_acetyl_p_pta;
    real coa_tilde = coa_sh_con            / km_coa_pta;
    
    
    real denom = (1 + acetyl_coa_tilde) * (1 + pi_tilde) + (1 + acetyl_p_tilde) * (1 + coa_tilde) - 1;
    
    real v = (kf_pta * acetyl_coa_tilde * pi_tilde - 
           kr_pta * acetyl_p_tilde * coa_tilde)/ denom * enzyme_pta;
    
    
    return v;
    
  }
  
  
  real compute_v_ack(real acetyl_p_con, real adp_con, real acetate_con, real atp_con, 
                            real km_acetyl_p_ack, real km_adp_ack, real km_acetate_ack, real km_atp_ack, 
                            real kf_ack, real kr_ack, real enzyme_ack){
  
    
    real acetyl_p_tilde = acetyl_p_con / km_acetyl_p_ack;
    real adp_tilde = adp_con / km_adp_ack;
    real atp_tilde = atp_con / km_atp_ack;
    real acetate_tilde = acetate_con / km_acetate_ack;
    
    
    real denom = (1 + acetyl_p_tilde) * (1 + adp_tilde) + (1 + acetate_tilde) * (1 + atp_tilde) - 1;
    
    real v = (kf_ack * acetyl_p_tilde * adp_tilde - kr_ack * acetate_tilde * atp_tilde) / denom * enzyme_ack;
    
    return v;
  
  }
  
  
    real compute_v_redox(real nadh_con, real nad_con, real k, real alpha){
  
  
    real v = k * alpha * nad_con - k * nadh_con;
    
    return v;
  
  }
  
  
  vector rhs(real t, vector y, array[] real theta, array[] real C, array[] real p_pdh,
      array[] real p_cs, array[] real p_aco, array[] real p_id, array[] real p_akd, array[] real p_scs,
      array[] real p_sd, array[] real p_fum, array[] real p_md, array[] real p_pk, array[] real p_pta,
      array[] real p_ack){
      
      
      vector[13] dy = rep_vector(0.0, 13);
  
      real NADH_c    = y[13];
      real NAD_c     = y[12];
      real ATP_c     = C[3];
      real ADP_c     = C[4];
      real FAD_c     = C[5];
      real FADH2_c   = C[6];
      real CO2_c     = C[7];
      real PEP_c     = C[8];
      real conserved_coa = C[9];
      real conserved_carbon = C[10];
      real pi_c = C[11];
      
      real CoA_free = conserved_coa - y[3] - y[8];
      real MAL = conserved_carbon - y[2] - y[4] - y[5] - y[6] - y[7] - y[8] - y[9] - y[10];



      real v_pdh = compute_v_pdh(
        y[1], CoA_free, NAD_c, y[3], NADH_c, CO2_c,
        // order: km_co2_pdh, km_pyruvate_pdh, km_coa_pdh, km_nad_pdh, km_acetyl_coa_pdh, km_nadh_pdh, enzyme_pdh, kf_pdh, kr_pdh, ki_nadh_pdh
        p_pdh[1], p_pdh[2], p_pdh[3] ,p_pdh[4], p_pdh[5], p_pdh[6], theta[1], p_pdh[7], p_pdh[9], p_pdh[8]);
  
      real v_cs = compute_v_cs(
        y[2], y[3], y[4], CoA_free, NADH_c,
        // order: km_oxalocetate_cs, km_acetyl_coa_cs, km_citrate_cs, km_coa_sh_cs, enzyme_cs, kf_cs, kr_cs, ki_nadh_cs
        p_cs[1], p_cs[2], p_cs[3], p_cs[4], theta[7], p_cs[5], p_cs[6], p_cs[7]);
  
      vector[2] v_aco = compute_v_aco(
        y[4], y[5], y[6],
        // order: km_citrate_aco, km_cis_aconitate_aco1, km_cis_aconitate_aco2,km_isocitrate_aco, enzyme_aco, kf1_aco, kr1_aco, kf2_aco, kr2_aco
        p_aco[1], p_aco[2], p_aco[3], p_aco[4], theta[2], p_aco[5], p_aco[6], p_aco[7], p_aco[8]);
      
      real v_id = compute_v_id(
        y[6], NAD_c, NADH_c, y[7], CO2_c,
        // order: km_co2_id, km_isocitrate_id, km_nadh_id, km_nad_id, km_aketo_id, enzyme_id, kf_id, kr_id
        p_id[1], p_id[2], p_id[3], p_id[4], p_id[5], theta[8], p_id[6], p_id[7]);
  
      real v_akd = compute_v_akd(
        y[7], NAD_c, CoA_free, y[8], NADH_c, CO2_c,
        // order: km_co2_akd, km_aketo_akd, km_nad_akd, km_nadh_akd, km_coa_sh_akd, km_succinyl_coa_akd, enzyme_akd, kf_akd, kr_akd
        p_akd[1], p_akd[2], p_akd[3], p_akd[4], p_akd[5], p_akd[6], theta[3], p_akd[7], p_akd[8]);
  
      real v_scs = compute_v_scs(
        y[8], ADP_c, y[9], CoA_free, NADH_c, ATP_c, pi_c,
        // order: km_succinyl_coa_scs, km_adp_scs, km_succinate_scs, km_coa_sh_scs, km_atp_scs, km_pi_scs, enzyme_scs, kf_scs, kr_scs
        p_scs[1], p_scs[2], p_scs[3], p_scs[4], p_scs[5], p_scs[8], theta[4], p_scs[6], p_scs[7]);
      
      real v_sd = compute_v_sd(
        y[9], FAD_c, y[10], FADH2_c,
        // order: km_succinate_sd, km_fad_sd, km_fumarate_sd, km_fadh2_sd, enzyme_sd, kf_sd, kr_sd
        p_sd[1], p_sd[2], p_sd[3], p_sd[4], theta[5], p_sd[5], p_sd[6]);
  
      real v_fum = compute_v_fumarase(
        y[10], MAL,
        // order: km_fum, km_mal, enz, kf, kr
        p_fum[1], p_fum[2], theta[12], p_fum[3], p_fum[4]);
  
      real v_md = compute_v_md(
        MAL, NAD_c, y[2], NADH_c,
        // order: km_malate_md, km_nad_md, km_oxalocetate_md, km_nadh_md, enzyme_md, kf_md, kr_md
        p_md[3], p_md[1], p_md[2], p_md[4], theta[11], p_md[5], p_md[6]);
  
      real v_pk = compute_v_pk(
      // km_pyruvate_pk, km_atp_pk, km_pep_pk, km_adp_pk, kf_pk, kr_pk, enzyme_pk, pyruvate, atp, pep, adp
        p_pk[5], p_pk[2], p_pk[1], p_pk[3], p_pk[4], p_pk[7], theta[6],
        y[1], ATP_c, PEP_c, ADP_c, p_pk[6]);
      
      
      real v_ack = compute_v_ack(
        //acetyl_p, adp, acetate, atp, km_acetyl_p_ack, km_adp_ack, km_acetate_ack, km_atp_ack, kf_ack, kr_ack, enzyme_ack
        y[11], ADP_c, theta[10], ATP_c, p_ack[1], p_ack[2], p_ack[3], p_ack[4], p_ack[5], p_ack[6], theta[13]);
    
    
      real v_pta = compute_v_pta(
        //acetyl_coa, pi, acetyl_p, coa_sh, km_acetyl_coa_pta, km_pi_pta, km_acetyl_p_pta, km_coa_pta, kf_pta, kr_pta, enzyme_pta
        y[3], pi_c, y[11], CoA_free, p_pta[1], p_pta[2], p_pta[3], p_pta[4], p_pta[5], p_pta[6], theta[9]);
    
  
    
      real v_redox = compute_v_redox(y[13], y[12], theta[14], theta[15]);
  
  
  
      dy[1]  +=  v_pk - v_pdh;            // PYR
      dy[2]  +=  v_md - v_cs;             // OAA
      dy[3]  +=  v_pdh - v_cs - v_pta;    // AcCoA
      dy[4]  +=  v_cs - v_aco[1];         // CIT
      dy[5]  +=  v_aco[1] - v_aco[2];     // CIS
      dy[6]  +=  v_aco[2] - v_id;         // ISO
      dy[7]  +=  v_id - v_akd;            // aKG
      dy[8]  +=  v_akd - v_scs;           // SucCoA
      dy[9]  +=  v_scs - v_sd;            // Succinate
      dy[10] +=  v_sd - v_fum;            // Fumarate
      dy[11] +=  v_pta - v_ack;            // Acetyl-P 
      dy[12] +=  - v_pdh - v_id - v_akd - v_md - v_redox; // NAD+
      dy[13] +=  v_pdh + v_id + v_akd + v_md + v_redox;   // NADH
   
      return dy;
  }
      
      
  array[] real compute_flux(vector y, array[] real theta, array[] real C, array[] real p_pdh,
    array[] real p_cs, array[] real p_aco, array[] real p_id, array[] real p_akd, array[] real p_scs,
    array[] real p_sd, array[] real p_fum, array[] real p_md, array[] real p_pk, array[] real p_pta,
    array[] real p_ack){
    
      array[14] real v;
      
      real NADH_c    = y[13];
      real NAD_c     = y[12];
      real ATP_c     = C[3];
      real ADP_c     = C[4];
      real FAD_c     = C[5];
      real FADH2_c   = C[6];
      real CO2_c     = C[7];
      real PEP_c     = C[8];
      real conserved_coa = C[9];
      real conserved_carbon = C[10];
      real pi_c = C[11];
      
      real CoA_free = conserved_coa - y[3] - y[8];
      real MAL = conserved_carbon - y[2] - y[4] - y[5] - y[6] - y[7] - y[8] - y[9] - y[10];

  
      real v_pdh = compute_v_pdh(
      y[1], CoA_free, NAD_c, y[3], NADH_c, CO2_c,
      // order: km_co2_pdh, km_pyruvate_pdh, km_coa_pdh, km_nad_pdh, km_acetyl_coa_pdh, km_nadh_pdh, enzyme_pdh, kf_pdh, kr_pdh, ki_nadh_pdh
      p_pdh[1], p_pdh[2], p_pdh[3] ,p_pdh[4], p_pdh[5], p_pdh[6], theta[1], p_pdh[7], p_pdh[9], p_pdh[8]);
  
      real v_cs = compute_v_cs(
      y[2], y[3], y[4], CoA_free, NADH_c,
      // order: km_oxalocetate_cs, km_acetyl_coa_cs, km_citrate_cs, km_coa_sh_cs, enzyme_cs, kf_cs, kr_cs, ki_nadh_cs
      p_cs[1], p_cs[2], p_cs[3], p_cs[4], theta[7], p_cs[5], p_cs[6], p_cs[7]);
  
      vector[2] v_aco = compute_v_aco(
      y[4], y[5], y[6],
      // order: km_citrate_aco, km_cis_aconitate_aco1, km_cis_aconitate_aco2,km_isocitrate_aco, enzyme_aco, kf1_aco, kr1_aco, kf2_aco, kr2_aco
      p_aco[1], p_aco[2], p_aco[3], p_aco[4], theta[2], p_aco[5], p_aco[6], p_aco[7], p_aco[8]);
      
      real v_id = compute_v_id(
      y[6], NAD_c, NADH_c, y[7], CO2_c,
      // order: km_co2_id, km_isocitrate_id, km_nadh_id, km_nad_id, km_aketo_id, enzyme_id, kf_id, kr_id
      p_id[1], p_id[2], p_id[3], p_id[4], p_id[5], theta[8], p_id[6], p_id[7]);
  
      real v_akd = compute_v_akd(
      y[7], NAD_c, CoA_free, y[8], NADH_c, CO2_c,
      // order: km_co2_akd, km_aketo_akd, km_nad_akd, km_nadh_akd, km_coa_sh_akd, km_succinyl_coa_akd, enzyme_akd, kf_akd, kr_akd
      p_akd[1], p_akd[2], p_akd[3], p_akd[4], p_akd[5], p_akd[6], theta[3], p_akd[7], p_akd[8]);
  
      real v_scs = compute_v_scs(
      y[8], ADP_c, y[9], CoA_free, NADH_c, ATP_c, pi_c,
      // order: km_succinyl_coa_scs, km_adp_scs, km_succinate_scs, km_coa_sh_scs, km_atp_scs, km_pi_scs, enzyme_scs, kf_scs, kr_scs
      p_scs[1], p_scs[2], p_scs[3], p_scs[4], p_scs[5], p_scs[8], theta[4], p_scs[6], p_scs[7]);
      
      real v_sd = compute_v_sd(
      y[9], FAD_c, y[10], FADH2_c,
      // order: km_succinate_sd, km_fad_sd, km_fumarate_sd, km_fadh2_sd, enzyme_sd, kf_sd, kr_sd
      p_sd[1], p_sd[2], p_sd[3], p_sd[4], theta[5], p_sd[5], p_sd[6]);
  
      real v_fum = compute_v_fumarase(
      y[10], MAL,
      // order: km_fum, km_mal, enz, kf, kr
      p_fum[1], p_fum[2], theta[12], p_fum[3], p_fum[4]);
  
      real v_md = compute_v_md(
      MAL, NAD_c, y[2], NADH_c,
      // order: km_malate_md, km_nad_md, km_oxalocetate_md, km_nadh_md, enzyme_md, kf_md, kr_md
      p_md[3], p_md[1], p_md[2], p_md[4], theta[11], p_md[5], p_md[6]);
  
      real v_pk = compute_v_pk(
      // km_pyruvate_pk, km_atp_pk, km_pep_pk, km_adp_pk, kf_pk, kr_pk, enzyme_pk, pyruvate, atp, pep, adp
      p_pk[5], p_pk[2], p_pk[1], p_pk[3], p_pk[4], p_pk[7], theta[6],
      y[1], ATP_c, PEP_c, ADP_c, p_pk[6]);
      
      
      real v_ack = compute_v_ack(
        //acetyl_p, adp, acetate, atp, km_acetyl_p_ack, km_adp_ack, km_acetate_ack, km_atp_ack, kf_ack, kr_ack, enzyme_ack
        y[11], ADP_c, theta[10], ATP_c, p_ack[1], p_ack[2], p_ack[3], p_ack[4], p_ack[5], p_ack[6], theta[13]);
    
    
      real v_pta = compute_v_pta(
        //acetyl_coa, pi, acetyl_p, coa_sh, km_acetyl_coa_pta, km_pi_pta, km_acetyl_p_pta, km_coa_pta, kf_pta, kr_pta, enzyme_pta
        y[3], pi_c, y[11], CoA_free, p_pta[1], p_pta[2], p_pta[3], p_pta[4], p_pta[5], p_pta[6], theta[9]);
  
  
      real v_redox = compute_v_redox(y[13], y[12], theta[14], theta[15]);
  
      // pack into output
      v[1]  = v_pk;
      v[2]  = v_pdh;
      v[3]  = v_cs;
      v[4]  = v_aco[1];
      v[5]  = v_aco[2];
      v[6]  = v_id;
      v[7]  = v_akd;
      v[8]  = v_scs;
      v[9]  = v_sd;
      v[10] = v_fum;
      v[11] = v_md;
      v[12] = v_pta;
      v[13] = v_ack;
      v[14] = v_redox;
  
      return v;
  }
      
      

}


data {
  int<lower=1> NY;
  int<lower=1> N_t;
  array[N_t] real t_obs;
  array[N_t, NY] real y_obs;
  real f_obs_ferm;
  real f_obs_tca;
  vector[NY] y0;

  array[11] real C;

  array[9] real p_pdh;
  array[7] real p_cs;
  array[8] real p_aco;
  array[7] real p_id;
  array[8] real p_akd;
  array[8] real p_scs;
  array[6] real p_sd;
  array[4] real p_fum;
  array[6] real p_md;
  array[7] real p_pk;
  array[6] real p_pta;
  array[6] real p_ack;
}
  
parameters {
    real pdh_raw;
    real aco_raw;
    real akd_raw;
    real scs_raw;
    real sd_raw;
    real pk_raw;
    real cs_raw;
    real id_raw;
    real pta_raw;
    real acetate_raw;
    real md_raw;
    real f_raw;
    real ack_raw;

    real k_raw;
    real alpha_raw;
    
    // raw noise parameters
    real sigma_c_raw;
    real sigma_f_raw;
}


transformed parameters {
    array[N_t] vector[NY] y_hat;
    array[14]       real f_hat;
    array[15]        real theta;
    real sigma_c;
    real sigma_f;
    
    // -----------------------------------------------
    theta[1] = pow(10, pdh_raw);   // PDH
    theta[2] = pow(10, aco_raw);   // ACO
    theta[3] = pow(10, akd_raw);   // AKGDH
    theta[4] = pow(10, scs_raw);   // SCS
    theta[5] = pow(10, sd_raw);    // SDH
    theta[6] = pow(10, pk_raw);    // PK 
    theta[7] = pow(10, cs_raw);    // CS
    theta[8] = pow(10, id_raw);    // IDH
    
    theta[9] = pow(10, pta_raw);    // CS
    theta[10] = 5;    // ACETATE
    theta[11] = pow(10, md_raw);    //MDH 0.023603447
    
    theta[12] = pow(10, f_raw); //0.001278269; //pow(10, f_raw);     // F
    theta[13] = pow(10, ack_raw); //0.002828854; //pow(10, ack_raw);   // ACKA
    
        
    theta[14] = pow(10, k_raw);         // k
    theta[15] = pow(10, alpha_raw);     // alpha
    
    sigma_c = pow(10, sigma_c_raw);  // 0.01–1
    sigma_f = pow(10, sigma_f_raw);  // pow(10, sigma_f_raw);

    // -----------------------------------------------
    // 4. ODE solve
    //    Relaxed tolerances to avoid CVode explosions
    // -----------------------------------------------
    y_hat = ode_bdf_tol(
    rhs, y0, 0.0, t_obs,
    1e-3, 1e-5, 1000000,
    theta, C, p_pdh, p_cs, p_aco, p_id, p_akd,
    p_scs, p_sd, p_fum, p_md, p_pk, p_pta, p_ack
    );

    f_hat = compute_flux(
    y_hat[N_t], theta, C, p_pdh, p_cs, p_aco, p_id, p_akd,
    p_scs, p_sd, p_fum, p_md, p_pk, p_pta, p_ack
    );
}

model {
    // priors — N(0,1) for all raw enzyme parameters
    pdh_raw ~ normal(-2,1);
    aco_raw ~ normal(-2,1); // normal(log10(0.005231645),0.3);
    akd_raw ~ normal(-2,1);
    scs_raw ~ normal(log10(0.0061704935), 0.5);
    sd_raw  ~ normal(-2,1);
    pk_raw  ~ normal(-2,1);
    cs_raw  ~ normal(-2,1);
    id_raw  ~ normal(-2,1);
    md_raw ~ normal(log10(0.023603447),0.5);
    f_raw ~ normal(log10(0.001278269),0.5);
    ack_raw ~ normal(-2,1);
    
    pta_raw  ~ normal(log10(0.000887043),0.5);
    acetate_raw  ~ normal(-1,1);
    
    k_raw ~ normal(log10(100),0.5);
    alpha_raw ~ normal(log10(0.0545/2.8),0.5);

    sigma_c_raw ~ normal(log10(0.2),0.15);
    sigma_f_raw ~ normal(log10(0.2),0.15);

    // likelihood: concentration
    for(j in 1:NY) {
        real yh = log10(fmax(y_hat[N_t][j], 1e-12));
        real yo = log10(fmax(y_obs[N_t,j], 1e-12));
        target += normal_lpdf(yo | yh, sigma_c);
    }


    // likelihood: tca flux
    for(j in 3:11){
        real fh_tca = log10(fmax(f_hat[j], 1e-12));
        real fo_tca = log10(fmax(f_obs_tca, 1e-12));
        
        target += normal_lpdf(fo_tca | fh_tca, sigma_f);
        
    }
    
    
    // likelihood: acetate flux
    for(j in 12:13){
        real fh_ferm = log10(fmax(f_hat[j], 1e-12));
        real fo_ferm = log10(fmax(f_obs_ferm, 1e-12));
        
        target += normal_lpdf(fo_ferm | fh_ferm, sigma_f);
    }
    
}




generated quantities {
  array[NY] real y_rep;
  array[14] real f_rep;

  for (j in 1:NY) {
    y_rep[j] = pow(
      10,
      normal_rng(log10(fmax(y_hat[1, j], 1e-12)), sigma_c)
    );
  }

  // Flux 可以为负，因此不要使用 log10 或 fmax
  for (j in 1:14) {
    f_rep[j] = normal_rng(f_hat[j], sigma_f);
  }
}



