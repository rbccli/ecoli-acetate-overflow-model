% clear;       
% clc;          
% close all; 

params = global_parameters_balanced();

% -------------------------------------------------------------------------
% Reference-baseline handling
%
% For the black (PDH + CS) reference scenario, set this to true. The script
% will draw the posterior samples, calculate a separate reference NADH and
% inhibition baseline for every sample, and save all paired quantities.
%
% For every other manually edited scenario, set this to false. The script
% will reload exactly the same theta rows and their matching baselines.
% -------------------------------------------------------------------------
compute_and_save_reference = false;
reference_file = 'posterior_reference_baselines_alpha.mat';

if compute_and_save_reference

    T = readtable('posterior1024.csv', ...
                  'VariableNamingRule', 'preserve');

    theta_all = table2array(T);

    % readtable can return a cell array when one or more CSV columns were
    % imported as text. Convert numeric-looking entries explicitly.
    if iscell(theta_all) || isstring(theta_all)
        theta_all = str2double(string(theta_all));
    end

    assert(isnumeric(theta_all) && ~iscell(theta_all), ...
        'Posterior samples could not be converted to a numeric matrix.');
    assert(all(isfinite(theta_all), 'all'), ...
        ['The posterior CSV contains missing or non-numeric entries. ', ...
         'Check the imported columns before running the analysis.']);
    assert(size(theta_all,2) >= 15, ...
        'The posterior file must contain at least the 15 expected columns.');

    rng(1);
    n_keep = 500;

    assert(size(theta_all,1) >= n_keep, ...
        'The posterior file contains fewer than %d rows.', n_keep);

    idx_keep = randperm(size(theta_all,1), n_keep);
    theta = theta_all(idx_keep, :);

    baseline_pdh_input = nan(n_keep,1);
    baseline_cs_input = nan(n_keep,1);
    nadh_reference_input = nan(n_keep,1);

else

    assert(isfile(reference_file), ...
        ['Reference file not found. First run the black PDH + CS scenario ', ...
         'with compute_and_save_reference = true.']);

    loaded = load(reference_file, 'reference');
    reference = loaded.reference;

    required_fields = {'theta','idx_keep','nadh_reference', ...
                       'baseline_pdh','baseline_cs'};

    assert(all(isfield(reference, required_fields)), ...
        'The saved reference file does not contain all required variables.');

    theta = reference.theta;
    idx_keep = reference.idx_keep;
    baseline_pdh_input = reference.baseline_pdh;
    baseline_cs_input = reference.baseline_cs;
    nadh_reference_input = reference.nadh_reference;

    n_keep = size(theta,1);

    assert(numel(baseline_pdh_input) == n_keep && ...
           numel(baseline_cs_input) == n_keep && ...
           numel(nadh_reference_input) == n_keep, ...
        'The saved theta samples and baseline vectors are not aligned.');

end


tspan = [0 60*60*24];


options = odeset('RelTol',1e-3,'AbsTol',1e-5);

y0 = [
    params.pyruvate_initial;
    params.oxaloacetate_initial;
    params.acetyl_coa_initial;
    params.citrate_initial;
    params.cis_aconitate_initial;
    params.isocitrate_initial;
    params.aketo_initial;
    params.succinyl_coa_initial;
    params.succinate_initial;
    params.fumarate_initial;
    params.acetyl_p_initial
    params.nad_pool * (1 - params.nadh_ratio);
    params.nad_pool * params.nadh_ratio
];

nadh_reference = nan(n_keep,1);
nadh_cs_baseline = nan(n_keep,1);
nadh_pdh_baseline = nan(n_keep,1);





fold_list = sort(unique([logspace(-1, log(20), 50), 1]));

n_grid = length(fold_list);

tca_flux  = nan(n_keep, n_grid);
ferm_flux = nan(n_keep, n_grid);
pdh_flux  = nan(n_keep, n_grid);

oaa_conc   = nan(n_keep, n_grid);
accoa_conc = nan(n_keep, n_grid);
pry_conc   = nan(n_keep, n_grid);
nad_conc   = nan(n_keep, n_grid);
nadh_conc  = nan(n_keep, n_grid);
alpha_list = nan(n_keep, n_grid);





parfor i = 1:n_keep


    p = params;

    p.enzyme_pdh = theta(i,1);
    p.enzyme_aco = theta(i,2);
    p.enzyme_akd = theta(i,3);
    p.enzyme_scs = theta(i,4);
    p.enzyme_sd  = theta(i,5);
    p.enzyme_pk  = theta(i,6);
    p.enzyme_cs  = theta(i,7);
    p.enzyme_id  = theta(i,8);
    p.enzyme_pta = theta(i,9);
    % p.acetate_constant = theta(i,10);
    p.enzyme_md  = theta(i,11);
    p.enzyme_fumarase   = theta(i,12);
    p.enzyme_ack = theta(i,13);

    p.k = theta(i,14);

    alpha = theta(i,15);   % posterior-fitted alpha
    p.alpha = alpha;       % reference state for this posterior sample

    if compute_and_save_reference

        % IMPORTANT: run this branch only while the rate functions contain
        % the full calibrated PDH + CS regulation (the black scenario).
        [~, reference_concentration] = ode15s( ...
            @(t, y) function_fsolve(y, p), ...
            tspan, y0, options);

        nadh_reference_i = reference_concentration(end,13);

        baseline_pdh_i = nadh_inhibition( ...
            nadh_reference_i, p.ki_nadh_pdh);

        baseline_cs_i = nadh_inhibition( ...
            nadh_reference_i, p.ki_nadh_cs);

    else

        % Reuse the baseline paired with this exact posterior sample.
        nadh_reference_i = nadh_reference_input(i);
        baseline_pdh_i = baseline_pdh_input(i);
        baseline_cs_i = baseline_cs_input(i);

    end

    nadh_reference(i) = nadh_reference_i;
    nadh_pdh_baseline(i) = baseline_pdh_i;
    nadh_cs_baseline(i) = baseline_cs_i;

    p.baseline_pdh = baseline_pdh_i;
    p.baseline_cs = baseline_cs_i;


    tca_row  = nan(1, n_grid);
    ferm_row = nan(1, n_grid);
    pdh_row  = nan(1, n_grid);

    oaa_row   = nan(1, n_grid);
    accoa_row = nan(1, n_grid);
    pry_row   = nan(1, n_grid);
    nad_row   = nan(1, n_grid);
    nadh_row  = nan(1, n_grid);

    for j = 1:n_grid

        p.alpha = fold_list(j) * alpha;

        % optional: define y0 consistently with alpha if alpha = NADH/NAD+
        y0_local = y0;

        [~, metabolite_concentration] = ode15s( ...
            @(t, y) function_fsolve(y, p), ...
            tspan, y0_local, options);

        fluxes = compute_flux(metabolite_concentration, p);

        if max(fluxes(end,3:11)) - min(fluxes(end,3:11)) <= 1e-4

            tca_row(j)  = fluxes(end,11);
            ferm_row(j) = fluxes(end,12);
            pdh_row(j)  = fluxes(end,1);

            oaa_row(j)   = metabolite_concentration(end,2);
            accoa_row(j) = metabolite_concentration(end,3);
            pry_row(j)   = metabolite_concentration(end,1);

            nad_row(j)  = metabolite_concentration(end,12);
            nadh_row(j) = metabolite_concentration(end,13);
        end

        % alpha_list(i,:) = fold_list(j) * k;
    end

    tca_flux(i,:)  = tca_row;
    ferm_flux(i,:) = ferm_row;
    pdh_flux(i,:)  = pdh_row;

    oaa_conc(i,:)   = oaa_row;
    accoa_conc(i,:) = accoa_row;
    pry_conc(i,:)   = pry_row;
    nad_conc(i,:)   = nad_row;
    nadh_conc(i,:)  = nadh_row;


    alpha_list(i,:) = fold_list * alpha;
end


%% Save the paired posterior samples and reference baselines
% This is executed only for the black PDH + CS reference scenario.
if compute_and_save_reference

    assert(all(isfinite(nadh_reference)), ...
        'At least one reference NADH value is non-finite.');
    assert(all(isfinite(nadh_pdh_baseline)) && ...
           all(nadh_pdh_baseline > 0), ...
        'At least one PDH baseline is invalid.');
    assert(all(isfinite(nadh_cs_baseline)) && ...
           all(nadh_cs_baseline > 0), ...
        'At least one CS baseline is invalid.');

    reference = struct();
    reference.theta = theta;
    reference.idx_keep = idx_keep(:);
    reference.nadh_reference = nadh_reference;
    reference.baseline_pdh = nadh_pdh_baseline;
    reference.baseline_cs = nadh_cs_baseline;
    reference.posterior_file = ...
        'posterior0914.csv';
    reference.rng_seed = 1;
    reference.n_keep = n_keep;
    reference.created_at = char(datetime('now', ...
        'Format', 'yyyy-MM-dd HH:mm:ss'));

    save(reference_file, 'reference', '-v7.3');

    fprintf(['Saved %d paired posterior samples and sample-specific ', ...
             'PDH/CS baselines to %s.\n'], ...
            n_keep, reference_file);

else

    fprintf(['Loaded %d paired posterior samples and sample-specific ', ...
             'PDH/CS baselines from %s.\n'], ...
            n_keep, reference_file);

end


%%
[fig1, ax1] = make_journal_figure(40,10);

% PDH + CS：黑色 
% PDH only：红色 -.
% CS only：蓝色
% No allostery：灰色 -.

black = [0   0   0]/255;        % PDH + CS
red   = [217 75 41]/255;        % PDH only
blue  = [8 127 163]/255;        % CS only
gray = [170 170 170]/255;       % No allostery

% figure(28);
ax2 = subplot(1,3,1);
plot_posterior_xy_common( ...
    alpha_list, ...
    tca_flux, ...
    black, ...
    '\alpha (Target NADH/NAD+)', ...
    'TCA cycle flux (mM s^{-1})', 'strict', 0.95, 50, 0.90);
% ylim([0.05 0.35]);
% yticks(linspace(0.1, 0.3, 3));
xlim([0.02 0.14]);
xticks([0.01 0.05 0.09 0.13]);

apply_journal_style(ax2, [], ...
    'AxisFontSize', 18, ...
    'LabelFontSize',18);



ax3 = subplot(1,3,2);
plot_posterior_xy_common( ...
    alpha_list, ...
    ferm_flux, ...
    black, ...
    '\alpha (Target NADH/NAD+)', ...
    'Pta-AckA flux (mM s^{-1})', 'strict', 0.95, 50, 0.90);
% ylim([0 1.25]);
% yticks(linspace(0.2, 1.2, 3));
xlim([0.02 0.14]);
xticks([0.01 0.05 0.09 0.13]);

apply_journal_style(ax3, [], ...
    'AxisFontSize', 18, ...
    'LabelFontSize', 18);



ratio = ferm_flux ./ (ferm_flux + tca_flux);

ax4 = subplot(1,3,3);
plot_posterior_xy_common( ...
    alpha_list, ...
    ratio, ...
    black, ...
    '\alpha (Target NADH/NAD+)', ...
    'Fermentation ratio', 'strict', 0.95, 50, 0.90);
% ylim([0 1]);
% yticks(linspace(0, 1, 3));
xlim([0.02 0.14]);
xticks([0.01 0.05 0.09 0.13]);

apply_journal_style(ax4, [], ...
    'AxisFontSize', 18, ...
    'LabelFontSize', 18);



%%

drawnow;

set([ax2, ax3, ax4], 'Units', 'normalized');

ax2.Position = [0.07, 0.20, 0.21, 0.72];
ax3.Position = [0.40, 0.20, 0.21, 0.72];
ax4.Position = [0.73, 0.20, 0.21, 0.72];

% export_figure(fig1, 'beta_inhibition_alpha0914');


%%



function plot_posterior_xy_common(X, Y, lineColor, xlab_txt, ylab_txt, rangeMode, coverage, n_grid, xCentralMass)

    if nargin < 6 || isempty(rangeMode)
        rangeMode = 'strict';
    end

    if nargin < 7 || isempty(coverage)
        coverage = 0.95;
    end

    if nargin < 8 || isempty(n_grid)
        n_grid = 50;
    end

    if nargin < 9 || isempty(xCentralMass)
        xCentralMass = 0.90;
    end

    valid = isfinite(X) & isfinite(Y);

    % Keep rows with at least two valid points
    valid_rows = sum(valid, 2) >= 2;

    X = X(valid_rows, :);
    Y = Y(valid_rows, :);
    valid = valid(valid_rows, :);

    % ============================================================
    % 1) Find x-range for each posterior sample
    % ============================================================

    n_samples = size(X, 1);

    row_x_min = nan(n_samples, 1);
    row_x_max = nan(n_samples, 1);

    for s = 1:n_samples

        xs = X(s, valid(s,:));

        if numel(xs) >= 2
            row_x_min(s) = min(xs);
            row_x_max(s) = max(xs);
        end

    end

    good_rows = isfinite(row_x_min) & ...
                isfinite(row_x_max) & ...
                row_x_max > row_x_min;

    X = X(good_rows, :);
    Y = Y(good_rows, :);
    valid = valid(good_rows, :);

    row_x_min = row_x_min(good_rows);
    row_x_max = row_x_max(good_rows);

    n_samples = size(X, 1);

    if n_samples < 2
        error('Not enough valid posterior samples to plot.');
    end

    % ============================================================
    % 2) Define common x-range
    % ============================================================

    switch lower(rangeMode)

        case 'strict'

            % Every posterior sample covers this range
            x_lo_common = max(row_x_min);
            x_hi_common = min(row_x_max);

        case 'robust'

            % Most posterior samples cover this range
            x_lo_common = prctile(row_x_min, coverage * 100);
            x_hi_common = prctile(row_x_max, (1 - coverage) * 100);

        otherwise

            error('rangeMode must be either ''strict'' or ''robust''.');

    end

    if x_hi_common <= x_lo_common
        error(['No common x-range across posterior samples. ', ...
               'Try rangeMode = ''robust'' or reduce coverage.']);
    end

    % ============================================================
    % 3) Further restrict to central xCentralMass of pooled x values
    % ============================================================

    x_all_valid = X(valid);
    x_all_valid = x_all_valid(isfinite(x_all_valid));

    if xCentralMass < 1

        tail = (1 - xCentralMass) / 2;

        x_lo_quantile = prctile(x_all_valid, 100 * tail);
        x_hi_quantile = prctile(x_all_valid, 100 * (1 - tail));

        % Final x-range = common support intersect central x range
        x_lo = max(x_lo_common, x_lo_quantile);
        x_hi = min(x_hi_common, x_hi_quantile);

    else

        x_lo = x_lo_common;
        x_hi = x_hi_common;

    end

    if x_hi <= x_lo
        error(['The intersection between common support and central x range is empty. ', ...
               'Try using xCentralMass = 1, or use rangeMode = ''robust''.']);
    end

    x_grid = linspace(x_lo, x_hi, n_grid);

    % ============================================================
    % 4) Interpolate every posterior trajectory onto common x-grid
    % ============================================================

    Y_grid = nan(n_samples, n_grid);

    for s = 1:n_samples

        xs = X(s, valid(s,:));
        ys = Y(s, valid(s,:));

        if numel(xs) < 2
            continue
        end

        % Sort by x
        [xs, ord] = sort(xs);
        ys = ys(ord);

        % Deal with duplicated x values by averaging y
        [xs_unique, ~, ic] = unique(xs);
        ys_unique = accumarray(ic(:), ys(:), [], @mean);

        if numel(xs_unique) < 2
            continue
        end

        % Only interpolate, never extrapolate
        in_range = x_grid >= min(xs_unique) & x_grid <= max(xs_unique);

        Y_grid(s, in_range) = interp1(xs_unique, ys_unique, ...
                                      x_grid(in_range), 'linear');

    end

    % ============================================================
    % 5) Compute median and quantile intervals at each x-grid point
    % ============================================================

    y_med = nan(1, n_grid);

    y_25 = nan(1, n_grid);
    y_75 = nan(1, n_grid);

    y_10 = nan(1, n_grid);
    y_90 = nan(1, n_grid);

    n_at_x = nan(1, n_grid);

    % Require enough samples at each x point
    min_required = max(8, ceil(0.8 * n_samples));

    for j = 1:n_grid

        yy = Y_grid(:, j);
        yy = yy(isfinite(yy));

        n_at_x(j) = numel(yy);

        if numel(yy) >= min_required

            y_med(j) = median(yy, 'omitnan');

            y_25(j) = prctile(yy, 25);
            y_75(j) = prctile(yy, 75);

            y_10(j) = prctile(yy, 10);
            y_90(j) = prctile(yy, 90);

        end

    end

    keep = isfinite(y_med);

    if ~any(keep)
        error('No x-grid points have enough posterior coverage.');
    end

    x_keep = x_grid(keep);

    y_med_keep = y_med(keep);

    y_25_keep = y_25(keep);
    y_75_keep = y_75(keep);

    y_10_keep = y_10(keep);
    y_90_keep = y_90(keep);

    % ============================================================
    % 6) Plot
    % ============================================================

    hold on

    fillColor = 1 - 0.45 * (1 - lineColor);

    % % 80% interval: 10% to 90%
    % fill([x_keep, fliplr(x_keep)], ...
    %      [y_10_keep, fliplr(y_90_keep)], ...
    %      fillColor, ...
    %      'EdgeColor', 'none', ...
    %      'FaceAlpha', 0.12);

    % 50% interval: 25% to 75%
    fill([x_keep, fliplr(x_keep)], ...
         [y_25_keep, fliplr(y_75_keep)], ...
         fillColor, ...
         'EdgeColor', 'none', ...
         'FaceAlpha', 0.12);

    % Plot a few posterior trajectories on the same x-grid
    % n_draw = min(5, n_samples);
    % idx_draw = randperm(n_samples, n_draw);
    % 
    % trajColor = 1 - 0.5 * (1 - lineColor);
    % 
    % for i = 1:n_draw
    % 
    %     s = idx_draw(i);
    % 
    %     ys = Y_grid(s, :);
    % 
    %     plot(x_grid, ys, '-', ...
    %         'Color', trajColor, ...
    %         'LineWidth', 0.8);
    % 
    % end

    % Median line
    plot(x_keep, y_med_keep, '-', ...
        'Color', lineColor, ...
        'LineWidth', 2.5);

    xlabel(xlab_txt);
    ylabel(ylab_txt);

    xlim([x_lo x_hi]);
    box on

end
