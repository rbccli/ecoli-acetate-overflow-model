function v = compute_v_pdh_beta(pyruvate_con, coa_con, nad_con, acetyl_coa_con, nadh_con, co2_con, km_co2, km_pyruvate, km_coa, km_nad, km_acetyl_coa, km_nadh, enzyme_pdh, kf_pdh, kr_pdh, ki_nadh, beta, baseline)
% Reaction 1: pyruvate + coa + nad -> co2 + h2o + acetyle-coa + nadh

    pyruvate_tilde = pyruvate_con / km_pyruvate;
    coa_tilde = coa_con / km_coa;
    nad_tilde = nad_con / km_nad;
    acetyl_coa_tilde = acetyl_coa_con/km_acetyl_coa;
    nadh_tilde = nadh_con / km_nadh;
    co2_tilde = co2_con/ km_co2;

    denominator = (1 + pyruvate_tilde) .* (1 + coa_tilde) .* (1 + nad_tilde) + (1 + acetyl_coa_tilde) .* (1 + nadh_tilde) .* (1 + co2_tilde) - 1;

    v = (kf_pdh .* pyruvate_tilde .* coa_tilde .* nad_tilde - kr_pdh .* acetyl_coa_tilde .* nadh_tilde .* co2_tilde)./ denominator .* enzyme_pdh .* baseline .* (nadh_inhibition(nadh_con, ki_nadh)./baseline).^ beta;

    % 0.0180112191884325; %

    % 0.6472; %
end