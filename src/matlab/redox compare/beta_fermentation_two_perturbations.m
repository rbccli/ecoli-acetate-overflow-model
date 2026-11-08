clear;

clc;

close all;





%% Load baseline parameters and posterior samples



params = global_parameters_balanced();



T = readtable( ...
    'posterior1024.csv');



theta = table2array(T);



rng(1);



n_keep = 1000;



assert(size(theta,1) >= n_keep, ...
    'n_keep exceeds the number of available posterior samples.');



assert(size(theta,2) >= 15, ...
    'The posterior file must contain at least 15 columns.');



idx_keep = randperm(size(theta,1), n_keep);
theta = theta(idx_keep, :);





%% ODE settings



tspan = [0, 60*60*24];



options = odeset('RelTol', 1e-3, ...
'AbsTol', 1e-5);





%% Initial conditions



% This assumes params.nadh_ratio is:

% NADH/(NADH + NAD+)

%

% If params.nadh_ratio represents NADH/NAD+, use:

%

% nadh_fraction = ...

%     params.nadh_ratio/(1 + params.nadh_ratio);



nadh_fraction = params.nadh_ratio;



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
    params.acetyl_p_initial;
    params.nad_pool * (1 - nadh_fraction);
    params.nad_pool * nadh_fraction
];





%% Perturbation grids



% Alpha fold change: alpha/alpha_0
alpha_fold_list = linspace(0.1, 10, 100);
alpha_fold_row  = alpha_fold_list(:)';



% NAD(H)-pool fold change: pool/pool_0
pool_fold_list = linspace(0.2, 5, 100);
pool_fold_row  = pool_fold_list(:)';



n_alpha_grid = numel(alpha_fold_list);
n_pool_grid  = numel(pool_fold_list);





%% Beta grid



beta_list = linspace(0, 1, 11);
n_beta = numel(beta_list);





%% Output arrays



% Fermentation slope versus alpha/alpha_0
alpha_ferm_beta_slope = nan(n_beta, 1);
alpha_ferm_slope_lo   = nan(n_beta, 1);
alpha_ferm_slope_hi   = nan(n_beta, 1);



% Fermentation slope versus NAD(H) pool/pool_0
pool_ferm_beta_slope = nan(n_beta, 1);
pool_ferm_slope_lo   = nan(n_beta, 1);
pool_ferm_slope_hi   = nan(n_beta, 1);



% Successful posterior-sample counts
n_valid_alpha = zeros(n_beta, 1);
n_valid_pool  = zeros(n_beta, 1);





%% Start parallel pool



if isempty(gcp('nocreate'))
    parpool;
end





%% Run beta scan with references anchored at beta = 1

fprintf('Inhibition references: beta = 1 for every posterior sample.\n');



for b = 1:n_beta



    beta_b = beta_list(b);



    fprintf( ...
        'Running beta = %.3f (%d/%d)\n', ...
        beta_b, b, n_beta);





    %% Store fermentation fluxes



    ferm_flux_alpha = nan(n_keep, n_alpha_grid);
    ferm_flux_pool  = nan(n_keep, n_pool_grid);





    parfor i = 1:n_keep



        %% Posterior-specific parameters



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
        p.acetate_constant = theta(i,10);
        p.enzyme_md  = theta(i,11);
        p.enzyme_fumarase   = theta(i,12);
        p.enzyme_ack = theta(i,13);
        p.k = theta(i,14);



        % Posterior baseline alpha

        alpha0 = theta(i,15);



        % Calculate this sample's inhibition references with full PDH feedback.
        % beta_b is applied only after the beta=1 references have been fixed.

        p.beta = 1;





        %% Sample-specific storage



        alpha_ferm_row = nan(1, n_alpha_grid);
        pool_ferm_row  = nan(1, n_pool_grid);





        %% Step 1: reference-state simulation at beta = 1



        p.alpha = alpha0;
        p.nad_pool = params.nad_pool;



        initial_baseline_successful = false;
        baseline_state_initial = nan(size(y0));



        try



            [~, baseline_solution_initial] = ode15s( ...
                @(t, y) function_fsolve(y, p), ...
                tspan, ...
                y0, ...
                options);



            baseline_state_initial = ...
                baseline_solution_initial(end,:)';



            baseline_dydt_initial = ...
                function_fsolve( ...
                baseline_state_initial, p);



            initial_baseline_successful = ...
                all(isfinite(baseline_state_initial)) && ...
                all(isfinite(baseline_dydt_initial)) && ...
                norm(baseline_dydt_initial, inf) <= 1e-6;



        catch

            initial_baseline_successful = false;

        end





        %% Step 2: fix inhibition references from the beta = 1 state



        if initial_baseline_successful



            baseline_nadh_initial = ...
                baseline_state_initial(13);



            p.baseline_pdh = nadh_inhibition( ...
                baseline_nadh_initial, ...
                p.ki_nadh_pdh);



            p.baseline_cs = nadh_inhibition( ...
                baseline_nadh_initial, ...
                p.ki_nadh_cs);



        end





        % Keep these beta=1 references throughout the current beta's baseline
        % and both perturbations. Identical beta=1 inputs give the same
        % sample-specific references at every beta in this scan.
        p.beta = beta_b;

        %% Step 3: verify the baseline at the scanned beta



        baseline_successful = false;
        baseline_state = nan(size(y0));



        if initial_baseline_successful && ...
                isfinite(p.baseline_pdh) && ...
                isfinite(p.baseline_cs) && ...
                p.baseline_pdh > 0 && ...
                p.baseline_cs > 0



            try



                [~, baseline_solution] = ode15s( ...
                    @(t, y) function_fsolve(y, p), ...
                    tspan, ...
                    baseline_state_initial, ...
                    options);



                baseline_state = ...
                    baseline_solution(end,:)';



                baseline_dydt = ...
                    function_fsolve(baseline_state, p);



                baseline_successful = ...
                    all(isfinite(baseline_state)) && ...
                    all(isfinite(baseline_dydt)) && ...
                    norm(baseline_dydt, inf) <= 1e-6;



            catch

                baseline_successful = false;

            end



        end





        %% Run both perturbation analyses



        if baseline_successful



            %% Baseline NAD(H) quantities



            baseline_nad  = baseline_state(12);

            baseline_nadh = baseline_state(13);



            baseline_pool = ...
                baseline_nad + baseline_nadh;



            valid_baseline_pool = ...
                isfinite(baseline_pool) && ...
                baseline_pool > 0;





            if valid_baseline_pool



                baseline_nad_fraction = ...
                    baseline_nad / baseline_pool;



                baseline_nadh_fraction = ...
                    baseline_nadh / baseline_pool;





                %% Perturbation 1: change alpha



                p_alpha = p;



                % Keep NAD(H) pool at baseline

                p_alpha.nad_pool = baseline_pool;



                for j = 1:n_alpha_grid



                    p_alpha.alpha = ...
                        alpha_fold_row(j) * alpha0;



                    % Begin at the baseline steady state

                    y0_alpha = baseline_state;



                    try



                        [~, alpha_solution] = ode15s( ...
                            @(t, y) ...
                            function_fsolve(y, p_alpha), ...
                            tspan, ...
                            y0_alpha, ...
                            options);



                        alpha_end = ...
                            alpha_solution(end,:)';



                        alpha_dydt = ...
                            function_fsolve( ...
                            alpha_end, p_alpha);



                        alpha_is_steady = ...
                            all(isfinite(alpha_end)) && ...
                            all(isfinite(alpha_dydt)) && ...
                            norm(alpha_dydt, inf) <= 1e-6;



                        if alpha_is_steady



                            alpha_fluxes = ...
                                compute_flux( ...
                                alpha_solution, p_alpha);



                            ferm_value = ...
                                alpha_fluxes(end,12);



                            if isfinite(ferm_value)



                                alpha_ferm_row(j) = ...
                                    ferm_value;



                            end



                        end



                    catch

                        % Failed simulation remains NaN

                    end



                end





                %% Perturbation 2: change NAD(H) pool size



                p_pool = p;



                % Keep alpha at its posterior baseline

                p_pool.alpha = alpha0;



                for j = 1:n_pool_grid



                    perturbed_pool = ...
                        pool_fold_row(j) * ...
                        baseline_pool;



                    % Update parameter-level pool

                    p_pool.nad_pool = perturbed_pool;



                    % Begin from baseline steady state

                    y0_pool = baseline_state;



                    % Change actual NAD+/NADH state-variable pool

                    % while preserving the baseline redox fraction

                    y0_pool(12) = ...
                        perturbed_pool * ...
                        baseline_nad_fraction;



                    y0_pool(13) = ...
                        perturbed_pool * ...
                        baseline_nadh_fraction;



                    try



                        [~, pool_solution] = ode15s( ...
                            @(t, y) ...
                            function_fsolve(y, p_pool), ...
                            tspan, ...
                            y0_pool, ...
                            options);



                        pool_end = ...
                            pool_solution(end,:)';



                        pool_dydt = ...
                            function_fsolve( ...
                            pool_end, p_pool);



                        pool_is_steady = ...
                            all(isfinite(pool_end)) && ...
                            all(isfinite(pool_dydt)) && ...
                            norm(pool_dydt, inf) <= 1e-6;



                        if pool_is_steady



                            pool_fluxes = ...
                                compute_flux( ...
                                pool_solution, p_pool);



                            ferm_value = ...
                                pool_fluxes(end,12);



                            if isfinite(ferm_value)



                                pool_ferm_row(j) = ...
                                    ferm_value;



                            end



                        end



                    catch

                        % Failed simulation remains NaN

                    end



                end



            end



        end





        %% Save this posterior sample



        ferm_flux_alpha(i,:) = alpha_ferm_row;

        ferm_flux_pool(i,:)  = pool_ferm_row;



    end





    %% Calculate one slope per posterior sample



    alpha_slopes_per_sample = nan(n_keep, 1);

    pool_slopes_per_sample  = nan(n_keep, 1);





    for i = 1:n_keep



        %% Fermentation slope versus alpha/alpha_0



        alpha_flux_i = ferm_flux_alpha(i,:)';



        valid_alpha_i = ...
            isfinite(alpha_fold_row(:)) & ...
            isfinite(alpha_flux_i);



        if sum(valid_alpha_i) >= 2



            alpha_coefficients = polyfit( ...
                alpha_fold_row(valid_alpha_i), ...
                alpha_flux_i(valid_alpha_i), ...
                1);



            alpha_slopes_per_sample(i) = ...
                alpha_coefficients(1);



        end





        %% Fermentation slope versus pool/pool_0



        pool_flux_i = ferm_flux_pool(i,:)';



        valid_pool_i = ...
            isfinite(pool_fold_row(:)) & ...
            isfinite(pool_flux_i);



        if sum(valid_pool_i) >= 2



            pool_coefficients = polyfit( ...
                pool_fold_row(valid_pool_i), ...
                pool_flux_i(valid_pool_i), ...
                1);



            pool_slopes_per_sample(i) = ...
                pool_coefficients(1);



        end



    end





    %% Filter unsuccessful samples



    valid_alpha_slopes = ...
        alpha_slopes_per_sample( ...
        isfinite(alpha_slopes_per_sample));



    valid_pool_slopes = ...
        pool_slopes_per_sample( ...
        isfinite(pool_slopes_per_sample));





    n_valid_alpha(b) = ...
        numel(valid_alpha_slopes);



    n_valid_pool(b) = ...
        numel(valid_pool_slopes);





    %% Posterior summaries



    if ~isempty(valid_alpha_slopes)



        alpha_ferm_beta_slope(b) = ...
            median(valid_alpha_slopes);



        alpha_ferm_slope_lo(b) = ...
            prctile(valid_alpha_slopes, 25);



        alpha_ferm_slope_hi(b) = ...
            prctile(valid_alpha_slopes, 75);



    end





    if ~isempty(valid_pool_slopes)



        pool_ferm_beta_slope(b) = ...
            median(valid_pool_slopes);



        pool_ferm_slope_lo(b) = ...
            prctile(valid_pool_slopes, 25);



        pool_ferm_slope_hi(b) = ...
            prctile(valid_pool_slopes, 75);



    end



end





%% Find zero-crossing beta values



beta_star_alpha = first_zero_crossing( ...
    beta_list, ...
    alpha_ferm_beta_slope);



beta_star_pool = first_zero_crossing( ...
    beta_list, ...
    pool_ferm_beta_slope);





fprintf('\nAlpha perturbation beta*: %.4f\n', ...
    beta_star_alpha);


fprintf('Pool perturbation beta*: %.4f\n', ...
    beta_star_pool);





%% Publication-style vertically aligned response-slope figure

[fig1, ax_dummy] = make_journal_figure(15, 14);

% Remove the single axis created by make_journal_figure
delete(ax_dummy);

tl = tiledlayout(fig1, 2, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax_alpha = nexttile(tl, 1);
ax_pool  = nexttile(tl, 2);


%% Colour palette

black = [0   0   0] / 255;
red   = [217 75 41] / 255;
blue  = [8 127 163] / 255;
gray  = [170 170 170] / 255;

feasible_colour      = [126 176 109] / 255;
feasible_text_colour = [75 110 65] / 255;

alpha_colour = red;
pool_colour  = blue;


%% Common beta-axis settings; y axes use independent ranges

common_xlim   = [-0.05, 1.05];
common_xticks = [0, 0.5, 1];


%% Convert plotting variables to row vectors

beta_plot = beta_list(:)';

alpha_slope_plot = alpha_ferm_beta_slope(:)';
alpha_lo_plot    = alpha_ferm_slope_lo(:)';
alpha_hi_plot    = alpha_ferm_slope_hi(:)';

pool_slope_plot = pool_ferm_beta_slope(:)';
pool_lo_plot    = pool_ferm_slope_lo(:)';
pool_hi_plot    = pool_ferm_slope_hi(:)';

% Include zero, the median, and the full uncertainty band in each panel.
alpha_ylim = independent_slope_ylim( ...
    [0, alpha_slope_plot, alpha_lo_plot, alpha_hi_plot]);
pool_ylim = independent_slope_ylim( ...
    [0, pool_slope_plot, pool_lo_plot, pool_hi_plot]);



%% Define thresholds and feasible region

have_alpha_star = ...
    isscalar(beta_star_alpha) && ...
    isfinite(beta_star_alpha);

have_pool_star = ...
    isscalar(beta_star_pool) && ...
    isfinite(beta_star_pool);

have_feasible_region = ...
    have_alpha_star && have_pool_star;

if have_feasible_region

    feasible_lower = min(beta_star_pool, beta_star_alpha);
    feasible_upper = max(beta_star_pool, beta_star_alpha);

end


%% ================================================================
%  Upper panel: target-redox perturbation
%  ================================================================

hold(ax_alpha, 'on');

ax_alpha.XScale = 'linear';

xlim(ax_alpha, common_xlim);
ylim(ax_alpha, alpha_ylim);

xticks(ax_alpha, common_xticks);
yticks(ax_alpha, 'auto');


%% Feasible region

if have_feasible_region

    h_feasible_patch = patch(ax_alpha, ...
        [feasible_lower, feasible_upper, ...
         feasible_upper, feasible_lower], ...
        [alpha_ylim(1), alpha_ylim(1), ...
         alpha_ylim(2), alpha_ylim(2)], ...
        feasible_colour, ...
        'FaceAlpha', 0.18, ...
        'EdgeColor', 'none');

end


%% Target-redox uncertainty interval

fill(ax_alpha, ...
    [beta_plot, fliplr(beta_plot)], ...
    [alpha_lo_plot, fliplr(alpha_hi_plot)], ...
    alpha_colour, ...
    'FaceAlpha', 0.14, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');


%% Target-redox median slope

plot(ax_alpha, ...
    beta_plot, ...
    alpha_slope_plot, ...
    '-', ...
    'Color', alpha_colour, ...
    'LineWidth', 2.2, ...
    'HandleVisibility', 'off');


%% Zero-response line

yline(ax_alpha, ...
    0, ...
    '--', ...
    'Color', gray, ...
    'LineWidth', 1.0, ...
    'HandleVisibility', 'off');


%% Threshold lines

% Pool threshold: lower boundary of the feasible region
if have_pool_star

    xline(ax_alpha, ...
        beta_star_pool, ...
        '--', ...
        'Color', pool_colour, ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

end


% Target-redox threshold with label
if have_alpha_star

    alpha_beta_label = sprintf( ...
        '\\beta^*_{\\alpha} = %.3f', ...
        beta_star_alpha);

    h_alpha_threshold = xline(ax_alpha, ...
        beta_star_alpha, ...
        '--', ...
        alpha_beta_label, ...
        'Color', alpha_colour, ...
        'LineWidth', 1.2, ...
        'FontSize', 12, ...
        'FontWeight', 'normal', ...
        'Interpreter', 'tex', ...
        'HandleVisibility', 'off');

    h_alpha_threshold.LabelVerticalAlignment = 'bottom';
    h_alpha_threshold.LabelHorizontalAlignment = 'right';

end


%% Upper-panel title

h_title_alpha = title(ax_alpha, ...
    'Target-redox perturbation', ...
    'FontWeight', 'normal');

% Only show beta tick labels on the lower panel
ax_alpha.XTickLabel = [];


%% ================================================================
%  Lower panel: NAD(H)-pool perturbation
%  ================================================================

hold(ax_pool, 'on');

ax_pool.XScale = 'linear';

xlim(ax_pool, common_xlim);
ylim(ax_pool, pool_ylim);

xticks(ax_pool, common_xticks);
yticks(ax_pool, 'auto');


%% Feasible region

if have_feasible_region

    patch(ax_pool, ...
        [feasible_lower, feasible_upper, ...
         feasible_upper, feasible_lower], ...
        [pool_ylim(1), pool_ylim(1), ...
         pool_ylim(2), pool_ylim(2)], ...
        feasible_colour, ...
        'FaceAlpha', 0.18, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

end


%% NAD(H)-pool uncertainty interval

fill(ax_pool, ...
    [beta_plot, fliplr(beta_plot)], ...
    [pool_lo_plot, fliplr(pool_hi_plot)], ...
    pool_colour, ...
    'FaceAlpha', 0.14, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');


%% NAD(H)-pool median slope

plot(ax_pool, ...
    beta_plot, ...
    pool_slope_plot, ...
    '-', ...
    'Color', pool_colour, ...
    'LineWidth', 2.2, ...
    'HandleVisibility', 'off');


%% Zero-response line

yline(ax_pool, ...
    0, ...
    '--', ...
    'Color', gray, ...
    'LineWidth', 1.0, ...
    'HandleVisibility', 'off');


%% Threshold lines

% Target-redox threshold: upper boundary of the feasible region
if have_alpha_star

    xline(ax_pool, ...
        beta_star_alpha, ...
        '--', ...
        'Color', alpha_colour, ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

end


% NAD(H)-pool threshold with label
if have_pool_star

    pool_beta_label = sprintf( ...
        '\\beta^*_{pool} = %.3f', ...
        beta_star_pool);

    h_pool_threshold = xline(ax_pool, ...
        beta_star_pool, ...
        '--', ...
        pool_beta_label, ...
        'Color', pool_colour, ...
        'LineWidth', 1.2, ...
        'FontSize', 12, ...
        'FontWeight', 'normal', ...
        'Interpreter', 'tex', ...
        'HandleVisibility', 'off');

    h_pool_threshold.LabelVerticalAlignment = 'bottom';
    h_pool_threshold.LabelHorizontalAlignment = 'left';

end


%% Lower-panel title and x-axis label

h_title_pool = title(ax_pool, ...
    'NAD(H)-pool perturbation', ...
    'FontWeight', 'normal');

h_xlabel = xlabel(ax_pool, ...
    'PDH inhibition scaling, \beta');


%% Apply journal style

apply_journal_style(ax_alpha, [], ...
    'AxisFontSize', 13, ...
    'LabelFontSize', 13);

apply_journal_style(ax_pool, [], ...
    'AxisFontSize', 13, ...
    'LabelFontSize', 13);


%% Restore axis settings after journal styling

for ax = [ax_alpha, ax_pool]

    ax.XScale = 'linear';

    xlim(ax, common_xlim);

    xticks(ax, common_xticks);

    ax.XColor = black;

    box(ax, 'off');

end

% Restore each panel's independent range after journal styling.
ylim(ax_alpha, alpha_ylim);
ylim(ax_pool, pool_ylim);
yticks(ax_alpha, 'auto');
yticks(ax_pool, 'auto');

% Suppress upper-panel x tick labels again
ax_alpha.XTickLabel = [];


%% Colour y axes according to perturbation type

ax_alpha.YColor = alpha_colour;
ax_pool.YColor  = pool_colour;


%% Shared quantity label; numerical y-axis ranges remain independent

h_ylabel_common = ylabel(tl, ...
    'Regression slope of fermentation flux');

set(h_ylabel_common, ...
    'FontName', ax_alpha.FontName, ...
    'FontSize', 13, ...
    'FontWeight', 'normal', ...
    'Color', black, ...
    'Interpreter', 'tex');


%% Explicit typography

set(h_xlabel, ...
    'FontName', ax_pool.FontName, ...
    'FontSize', 13, ...
    'FontWeight', 'normal', ...
    'Interpreter', 'tex');

set([h_title_alpha, h_title_pool], ...
    'FontName', ax_alpha.FontName, ...
    'FontSize', 13, ...
    'FontWeight', 'normal', ...
    'Color', black, ...
    'Interpreter', 'tex');


%% Preserve threshold-label typography

if exist('h_alpha_threshold', 'var') && ...
        isgraphics(h_alpha_threshold)

    h_alpha_threshold.FontName   = ax_alpha.FontName;
    h_alpha_threshold.FontSize   = 12;
    h_alpha_threshold.FontWeight = 'normal';

end

if exist('h_pool_threshold', 'var') && ...
        isgraphics(h_pool_threshold)

    h_pool_threshold.FontName   = ax_pool.FontName;
    h_pool_threshold.FontSize   = 12;
    h_pool_threshold.FontWeight = 'normal';

end


%% Legend explaining the common feasible region

if have_feasible_region

    h_feasible_legend = legend( ...
        ax_alpha, ...
        h_feasible_patch, ...
        'Feasible for both responses', ...
        'Location', 'northeast', ...
        'Box', 'off', ...
        'FontSize', 10, ...
        'FontWeight', 'normal', ...
        'Interpreter', 'tex');

    h_feasible_legend.TextColor = feasible_text_colour;

end


%% Link the beta axes

linkaxes([ax_alpha, ax_pool], 'x');

xlim(ax_alpha, common_xlim);
xlim(ax_pool, common_xlim);

drawnow;

%% Export

% fig1.Color = 'white';
% 
% exportgraphics( ...
%     fig1, ...
%     'beta_slope0914.pdf', ...
%     'ContentType', 'vector', ...
%     'BackgroundColor', 'white');

% exportgraphics( ...
%     fig1, ...
%     'beta_slope.png', ...
%     'Resolution', 600, ...
%     'BackgroundColor', 'white');



%% Export figure



% export_figure( ...
%     fig1, ...
%     'beta_fermentation_two_perturbations');





%% Save numerical results



results_table = table( ...
    beta_list', ...
    alpha_ferm_beta_slope, ...
    alpha_ferm_slope_lo, ...
    alpha_ferm_slope_hi, ...
    pool_ferm_beta_slope, ...
    pool_ferm_slope_lo, ...
    pool_ferm_slope_hi, ...
    n_valid_alpha, ...
    n_valid_pool, ...
    'VariableNames', { ...
        'Beta', ...
        'AlphaSlopeMedian', ...
        'AlphaSlopeLower50', ...
        'AlphaSlopeUpper50', ...
        'PoolSlopeMedian', ...
        'PoolSlopeLower50', ...
        'PoolSlopeUpper50', ...
        'NValidAlpha', ...
        'NValidPool'});


writetable( ...
    results_table, ...
    'beta_fermentation_two_perturbations.csv');





%% Local function



function beta_star = first_zero_crossing( ...
    beta_values, slope_values)



    beta_values = beta_values(:);

    slope_values = slope_values(:);



    valid = ...
        isfinite(beta_values) & ...
        isfinite(slope_values);



    x = beta_values(valid);

    y = slope_values(valid);



    beta_star = NaN;



    if numel(x) < 2

        return;

    end



    exact_zero = find(y == 0, 1, 'first');



    if ~isempty(exact_zero)



        beta_star = x(exact_zero);

        return;



    end



    crossing_index = find( ...
        y(1:end-1) .* y(2:end) < 0, ...
        1, ...
        'first');



    if isempty(crossing_index)

        return;

    end



    x1 = x(crossing_index);

    x2 = x(crossing_index + 1);



    y1 = y(crossing_index);

    y2 = y(crossing_index + 1);



    beta_star = ...
        x1 - y1 * (x2 - x1)/(y2 - y1);



end


function limits = independent_slope_ylim(values)
    values = values(isfinite(values));
    if isempty(values)
        limits = [-0.01, 0.01];
        return;
    end
    lower = min(values(:));
    upper = max(values(:));
    padding = max(0.12 * (upper - lower), 1e-5);
    limits = [lower - padding, upper + padding];
end
