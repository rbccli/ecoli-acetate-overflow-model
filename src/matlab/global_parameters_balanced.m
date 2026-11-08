function params = global_parameters_balanced()

    params = struct();


    % inihibition constant
    
    params.ki_nadh_pdh = 0.001;
    params.ki_nadh_cs = 0.0028;
    params.ki_atp_pk = 4;


    % michaelis constant
    % pk
    params.km_pyruvate_pk = 0.0505;
    params.km_pep_pk = 0.2666;
    params.km_atp_pk = 16.1752;
    params.km_adp_pk = 0.1919;


    % pdh
    params.km_pyruvate_pdh = 0.3296; 
    params.km_coa_pdh = 0.0471; 
    params.km_nad_pdh = 0.0619; 
    params.km_nadh_pdh = 0.215;
    params.km_acetyl_coa_pdh = 0.215;
    params.km_co2_pdh = 0.215;

    % pta
    params.km_acetyl_coa_pta = 0.0491;
    params.km_pi_pta = 1.9606;
    params.km_acetyl_p_pta = 0.6297;
    params.km_coa_pta = 0.0521;
    
    % ack
    params.km_acetyl_p_ack = 0.0938;
    params.km_adp_ack = 0.2799;
    params.km_acetate_ack = 7.9678;
    params.km_atp_ack = 0.096;
    

    % cs
    params.km_citrate_cs = 0.0018; 
    params.km_coa_sh_cs = 44914.6131; % 0.0027; 
    params.km_oxalocetate_cs = 0.0119; 
    params.km_acetyl_coa_cs = 0.0521; 


    % aconitase
    params.km_citrate_aco = 0.6464; 
    params.km_cis_aconitate_aco1 = 0.0095;
    params.km_cis_aconitate_aco2 = 0.2098;
    params.km_isocitrate_aco = 0.0872; 

    % isocitrate dehydrogenase
    params.km_isocitrate_id = 0.0109; 
    params.km_nad_id = 1.5211; 
    params.km_aketo_id = 0.4744; 
    params.km_nadh_id = 0.007; % 24.3; 
    params.km_co2_id = 2.4317;


    % akd checked not much different
    params.km_nad_akd = 0.0512; 
    params.km_aketo_akd = 0.0026; 
    params.km_coa_sh_akd = 0.0512; 
    params.km_nadh_akd = 0.1953; 
    params.km_succinyl_coa_akd = 0.1953; 
    params.km_co2_akd = 0.1953;


    % scs finished check 
    params.km_coa_sh_scs = 2.0659; 
    params.km_succinyl_coa_scs = 0.0316; 
    params.km_adp_scs = 0.3244; 
    params.km_succinate_scs = 0.1444; 
    params.km_atp_scs = 0.0426; 
    params.km_pi_scs = 0.8953;


    % sd finished check 
    params.km_succinate_sd = 0.0021;  
    params.km_q_sd = 0.002; % q
    params.km_fumarate_sd = 0.0233; 
    params.km_qh2_sd = 12.5427; % qh2


    % fumarase
    params.km_fumarate_famarase = 0.4222; 
    params.km_malate_fumarase = 0.6905; 

    % malate dehydrogenase dont change same ecoli exp
    params.km_nad_md = 0.418; 
    params.km_oxalocetate_md = 0.0244; 
    params.km_malate_md = 3.8086; 
    params.km_nadh_md = 0.0301; 





    % params.enzyme_pk = 0.287; %0.460;
    % params.enzyme_pdh = 2.39; %0.564;
    % params.enzyme_cs = 3.81; %0.514;
    % params.enzyme_aco = 2.14;
    % params.enzyme_id = 0.0268; %0.000536;
    % params.enzyme_akd = 0.167; %0.00202;
    % params.enzyme_scs = 0.262; %0.0476;
    % params.enzyme_sd = 0.00907; %0.000177;
    % params.enzyme_fumarase = 0.00325;
    % params.enzyme_md = 0.0471;
    % params.enzyme_pta = 0.488; %0.00140;
    % params.enzyme_ack = 0.00197;



    params.acetate_constant = 5;  %4.49;



    params.enzyme_pdh = 2.11;
    params.enzyme_aco = 1.21;
    params.enzyme_akd = 0.208;
    params.enzyme_scs = 0.0277;
    params.enzyme_sd  = 0.00911;
    params.enzyme_pk  = 0.0785;
    params.enzyme_cs  = 15.5;
    params.enzyme_id  = 0.0190;
    params.enzyme_pta = 0.0589;
    params.enzyme_md  = 0.0534;
    params.enzyme_fumarase   = 0.00412;
    params.enzyme_ack = 0.00513;




    % forward and reverse k_cat 
    params.kf_pk = 555.4713;
    params.kr_pk = 0.4124;
    params.kf_pdh = 36.1192; 
    params.kr_pdh = 0;
    params.kf_cs = 141.5176; 
    params.kr_cs = 5.8632; 
    params.kf1_aco = 5.7836; 
    params.kf2_aco = 4.5447; 
    params.kr1_aco = 3.8424; 
    params.kr2_aco = 1.0233; 
    params.kf_id = 82.0188; 
    params.kr_id = 14.7374;
    params.kf_akd = 2.0345; 
    params.kr_akd = 0.0011;
    params.kf_scs = 35.8092; 
    params.kr_scs = 60.5031; 
    params.kf_sd = 32.8092; 
    params.kr_sd = 180.5162; 
    params.kf_fumarase = 723.0053; 
    params.kr_fumarase = 395.3746;
    params.kf_md = 132.9511; 
    params.kr_md = 1557.9368; 
    params.kf_pta = 22.1455; 
    params.kr_pta = 237.7203;
    params.kf_ack = 4791.0944;
    params.kr_ack = 852.687;


    params.pyruvate_initial = 3.66;
    params.oxaloacetate_initial = 0.000487;
    params.acetyl_coa_initial = 0.606;
    params.citrate_initial = 1.96;
    params.cis_aconitate_initial = 0.0161;
    params.isocitrate_initial = 0.0367;
    params.aketo_initial = 0.443;
    params.coa_initial = 1.37;
    params.succinyl_coa_initial = 0.233;
    params.succinate_initial = 0.569;
    params.fumarate_initial = 0.115;
    params.malate_initial = 1.68;

    params.acetyl_p_initial = 1.07;


    
    % params.nadh_constant = 0.0545; %% lower bound value
    % params.atp_constant = 8.13;
    % params.qh2_constant = 0.0173; %% qh2 no data 
    % 
    % 
    % params.adp_constant = 0.704;
    % params.q_constant = 0.173; %% q no data
    % params.nad_constant = 2.80;


    params.co2_constant = 0.0752; 
    params.pi_constant = 23.9;
    params.pep_constant = 0.231; %0.184;


    params.adp_pool = 8.8340; % params.adp_constant + params.atp_constant;
    params.q_pool = 0.1903; % params.qh2_constant + params.q_constant;
    params.nad_pool = 2.8545; % params.nad_constant + params.nadh_constant;

    params.atp_ratio = 0.9203; %params.atp_constant / params.adp_pool;
    params.qh2_ratio = 0.0909; %params.qh2_constant / params.q_pool;
    params.nadh_ratio = 0.0191; %params.nadh_constant / params.nad_pool;



    params.conserved_coa = params.coa_initial + params.succinyl_coa_initial + params.acetyl_coa_initial;
    params.conserved_quantity2 = params.oxaloacetate_initial + params.citrate_initial + params.cis_aconitate_initial + params.isocitrate_initial + params.aketo_initial + params.succinyl_coa_initial + params.succinate_initial + params.fumarate_initial + params.malate_initial;



    params.k = 187; %10000;
    params.alpha = 0.00627; %0.0545/2.80;


    params.beta = 0.08;
    params.baseline_pdh = 0.0180112191884325;
    params.baseline_cs = 1;
end



