%% Sobol sensitivity analysis using total-order indices
clear;
clc;
close all;

rng default;   % For reproducibility


%% ================================================================
%  1. Load baseline parameters and choose parameters to vary
%  ================================================================

params = global_parameters_balanced();

fields = fieldnames(params);

fields = fields(~( ...
    contains(fields, 'initial', 'IgnoreCase', true) | ...
    contains(fields, 'nadh_ratio', 'IgnoreCase', true) | ...
    contains(fields, 'atp', 'IgnoreCase', true) | ...
    contains(fields, 'adp', 'IgnoreCase', true) | ...
    contains(fields, 'acetate_constant', 'IgnoreCase', true) | ...
    contains(fields, 'co2_constant', 'IgnoreCase', true) | ...
    contains(fields, 'pi_constant', 'IgnoreCase', true) | ...
    contains(fields, 'pep_constant', 'IgnoreCase', true) | ...
    contains(fields, 'q', 'IgnoreCase', true) | ...
    contains(fields, 'km', 'IgnoreCase', true) | ...
    contains(fields, 'conserved', 'IgnoreCase', true) | ...
    contains(fields, 'beta', 'IgnoreCase', true) | ...
    contains(fields, 'kf', 'IgnoreCase', true) | ...
    contains(fields, 'baseline', 'IgnoreCase', true) | ...
    contains(fields, 'kr', 'IgnoreCase', true) ...
));

num_samples = 10000;
num_params  = numel(fields);

num_metabolite_outputs = 13;
num_total_outputs      = 15;   % 13 states + 2 fluxes


%% Extract baseline parameter values

baseline = zeros(1, num_params);

for j = 1:num_params
    baseline(j) = params.(fields{j});
end

if any(~isfinite(baseline)) || any(baseline <= 0)
    bad_parameters = fields(~isfinite(baseline) | baseline <= 0);

    error( ...
        ['All parameters must be finite and positive for log-space ', ...
         'sampling. Problematic parameters:\n%s'], ...
        strjoin(bad_parameters, ', ') ...
    );
end

fprintf('Number of varied parameters: %d\n', num_params);
fprintf('Requested Sobol base samples: %d\n', num_samples);


%% ================================================================
%  2. Generate A and B matrices in log10 space
%  ================================================================

% Each parameter varies from 0.1 times to 10 times its baseline value
%% Generate A and B matrices exactly as in the original code

Varmin = baseline * 0.1;
Varmax = baseline * 10;

Varmin_AB = [Varmin, Varmin];
Varmax_AB = [Varmax, Varmax];

% This reproduces the Sobol sequence actually used in the original code
p = sobolset(2 * num_params);

U = p(1:num_samples, :);

logVarmin = log10(Varmin_AB);
logVarmax = log10(Varmax_AB);

log_R = logVarmin + U .* (logVarmax - logVarmin);

R = 10.^log_R;

samples_A = R(:, 1:num_params);
samples_B = R(:, num_params + 1:end);

fprintf('Sobol sampling completed.\n');


%% ================================================================
%  3. Construct parameter structures for A and B
%  ================================================================

Aparams = repmat(params, num_samples, 1);
Bparams = repmat(params, num_samples, 1);

for j = 1:num_params

    fname = fields{j};

    values_A = num2cell(samples_A(:,j));
    [Aparams.(fname)] = values_A{:};

    values_B = num2cell(samples_B(:,j));
    [Bparams.(fname)] = values_B{:};

end


%% Recalculate conserved pools

for i = 1:num_samples

    Aparams(i).conserved_coa = ...
        Aparams(i).coa_initial + ...
        Aparams(i).succinyl_coa_initial + ...
        Aparams(i).acetyl_coa_initial;

    Aparams(i).conserved_quantity2 = ...
        Aparams(i).oxaloacetate_initial + ...
        Aparams(i).citrate_initial + ...
        Aparams(i).cis_aconitate_initial + ...
        Aparams(i).isocitrate_initial + ...
        Aparams(i).aketo_initial + ...
        Aparams(i).succinyl_coa_initial + ...
        Aparams(i).succinate_initial + ...
        Aparams(i).fumarate_initial + ...
        Aparams(i).malate_initial;


    Bparams(i).conserved_coa = ...
        Bparams(i).coa_initial + ...
        Bparams(i).succinyl_coa_initial + ...
        Bparams(i).acetyl_coa_initial;

    Bparams(i).conserved_quantity2 = ...
        Bparams(i).oxaloacetate_initial + ...
        Bparams(i).citrate_initial + ...
        Bparams(i).cis_aconitate_initial + ...
        Bparams(i).isocitrate_initial + ...
        Bparams(i).aketo_initial + ...
        Bparams(i).succinyl_coa_initial + ...
        Bparams(i).succinate_initial + ...
        Bparams(i).fumarate_initial + ...
        Bparams(i).malate_initial;

end

fprintf('Parameter structures completed.\n');


%% ================================================================
%  4. Start parallel pool
%  ================================================================

pool = gcp('nocreate');

if isempty(pool)

    parpool('threads');

elseif ~isa(pool, 'parallel.ThreadPool')

    delete(pool);
    parpool('threads');

end


%% ================================================================
%  5. Simulate A, B and all AB hybrid matrices
%  ================================================================

%% ================================================================
%  5. Simulate in batches with checkpointing
%  ================================================================

steady_YA = nan(num_samples, num_metabolite_outputs);
steady_YB = nan(num_samples, num_metabolite_outputs);

steady_YAB = nan( ...
    num_samples, ...
    num_params, ...
    num_metabolite_outputs ...
);

okA  = false(num_samples, 1);
okB  = false(num_samples, 1);
okAB = false(num_samples, num_params);

% Records completed outer parfor iterations, including failed ODE runs
completed_rows = false(num_samples, 1);

batch_size = 500;

checkpoint_file = 'sobol_checkpoint.mat';
checkpoint_tmp  = 'sobol_checkpoint_tmp.mat';

checkpoint_num_samples = num_samples;
checkpoint_fields       = fields;
checkpoint_baseline     = baseline;


%% Evaluate all samples from scratch
% Saved checkpoints are outputs only and are never loaded on startup.

sample_indices = (1:num_samples)';
num_batches = ceil(num_samples / batch_size);

fprintf( ...
    'Beginning simulations from scratch: %d samples in %d batches.\n', ...
    num_samples, ...
    num_batches ...
);

session_timer = tic;

%% Run batches

for batch_number = 1:num_batches

    first_position = (batch_number - 1) * batch_size + 1;

    last_position = min( ...
        batch_number * batch_size, ...
        num_samples ...
    );

    batch_indices = sample_indices( ...
        first_position:last_position ...
    );

    current_batch_size = numel(batch_indices);


    % Batch-local output arrays

    batch_YA = nan( ...
        current_batch_size, ...
        num_metabolite_outputs ...
    );

    batch_YB = nan( ...
        current_batch_size, ...
        num_metabolite_outputs ...
    );

    batch_YAB = nan( ...
        current_batch_size, ...
        num_params, ...
        num_metabolite_outputs ...
    );

    batch_okA  = false(current_batch_size, 1);
    batch_okB  = false(current_batch_size, 1);
    batch_okAB = false(current_batch_size, num_params);


    batch_timer = tic;


    parfor kk = 1:current_batch_size

        i = batch_indices(kk);

        local_YAB  = nan(num_params, num_metabolite_outputs);
        local_okAB = false(1, num_params);


        % A simulation

        [local_okA, yAend] = run_one(Aparams(i));

        batch_okA(kk)  = local_okA;
        batch_YA(kk,:) = yAend;


        % B and AB simulations only when required

        if local_okA

            [local_okB, yBend] = run_one(Bparams(i));

            batch_okB(kk)  = local_okB;
            batch_YB(kk,:) = yBend;

            if local_okB

                for h = 1:num_params

                    pAB = Aparams(i);

                    fname = fields{h};
                    pAB.(fname) = Bparams(i).(fname);

                    [local_okAB(h), yABend] = run_one(pAB);

                    local_YAB(h,:) = yABend;

                end

            end

        end

        batch_okAB(kk,:) = local_okAB;

        batch_YAB(kk,:,:) = reshape( ...
            local_YAB, ...
            1, ...
            num_params, ...
            num_metabolite_outputs ...
        );

    end


    %% Copy completed batch into the full arrays

    steady_YA(batch_indices,:) = batch_YA;
    steady_YB(batch_indices,:) = batch_YB;
    steady_YAB(batch_indices,:,:) = batch_YAB;

    okA(batch_indices) = batch_okA;
    okB(batch_indices) = batch_okB;
    okAB(batch_indices,:) = batch_okAB;

    % Mark completed even when an individual ODE failed
    completed_rows(batch_indices) = true;


    %% Save checkpoint safely

    save( ...
        checkpoint_tmp, ...
        'steady_YA', ...
        'steady_YB', ...
        'steady_YAB', ...
        'okA', ...
        'okB', ...
        'okAB', ...
        'completed_rows', ...
        'checkpoint_num_samples', ...
        'checkpoint_fields', ...
        'checkpoint_baseline', ...
        '-v7.3' ...
    );

    [move_success, move_message] = movefile( ...
        checkpoint_tmp, ...
        checkpoint_file, ...
        'f' ...
    );

    if ~move_success
        error('Could not update checkpoint: %s', move_message);
    end


    %% Report progress and estimated remaining time

    total_completed = nnz(completed_rows);

    elapsed_minutes = toc(session_timer) / 60;

    samples_per_minute = ...
        total_completed / max(elapsed_minutes, eps);

    estimated_minutes_remaining = ...
        (num_samples - total_completed) / samples_per_minute;

    fprintf( ...
        ['Batch %d/%d completed in %.2f min | ', ...
         '%d/%d samples (%.1f%%) | ETA %.1f min\n'], ...
        batch_number, ...
        num_batches, ...
        toc(batch_timer) / 60, ...
        total_completed, ...
        num_samples, ...
        100 * total_completed / num_samples, ...
        estimated_minutes_remaining ...
    );

end


simulation_time = toc(session_timer);

fprintf( ...
    'Simulation stage completed in %.1f minutes this session.\n', ...
    simulation_time / 60 ...
);


%% ================================================================
%  6. First filter: successful ODE integrations
%  ================================================================

ode_good = okA & okB & all(okAB, 2);

IA = find(ode_good);

fprintf( ...
    'Samples passing ODE integration: %d/%d\n', ...
    numel(IA), ...
    num_samples ...
);

if isempty(IA)
    error('No sample passed the ODE integration criteria.');
end

steady_YA  = steady_YA(IA,:);
steady_YB  = steady_YB(IA,:);
steady_YAB = steady_YAB(IA,:,:);

Aparams = Aparams(IA);
Bparams = Bparams(IA);

num_ode_valid = numel(IA);


%% ================================================================
%  7. Calculate fluxes and apply flux constraints
%  ================================================================

tca_flux_A  = nan(num_ode_valid, 1);
tca_flux_B  = nan(num_ode_valid, 1);
tca_flux_AB = nan(num_ode_valid, num_params);

ferm_flux_A  = nan(num_ode_valid, 1);
ferm_flux_B  = nan(num_ode_valid, 1);
ferm_flux_AB = nan(num_ode_valid, num_params);

valid_flux_A  = false(num_ode_valid, 1);
valid_flux_B  = false(num_ode_valid, 1);
valid_flux_AB = false(num_ode_valid, num_params);

fprintf('Evaluating flux constraints...\n');

parfor i = 1:num_ode_valid

    local_tca_AB   = nan(1, num_params);
    local_ferm_AB  = nan(1, num_params);
    local_valid_AB = false(1, num_params);


    % --------------------------------------------------------------
    % A matrix
    % --------------------------------------------------------------

    yA = steady_YA(i,:);

    [ ...
        valid_flux_A(i), ...
        tca_flux_A(i), ...
        ferm_flux_A(i) ...
    ] = extract_flux_outputs(yA, Aparams(i));


    % --------------------------------------------------------------
    % B matrix
    % --------------------------------------------------------------

    yB = steady_YB(i,:);

    [ ...
        valid_flux_B(i), ...
        tca_flux_B(i), ...
        ferm_flux_B(i) ...
    ] = extract_flux_outputs(yB, Bparams(i));


    % --------------------------------------------------------------
    % AB hybrid matrices
    % --------------------------------------------------------------

    for h = 1:num_params

        pAB = Aparams(i);

        fname = fields{h};

        pAB.(fname) = Bparams(i).(fname);

        yAB = reshape( ...
            steady_YAB(i,h,:), ...
            1, ...
            num_metabolite_outputs ...
        );

        [ ...
            local_valid_AB(h), ...
            local_tca_AB(h), ...
            local_ferm_AB(h) ...
        ] = extract_flux_outputs(yAB, pAB);

    end

    valid_flux_AB(i,:) = local_valid_AB;
    tca_flux_AB(i,:)   = local_tca_AB;
    ferm_flux_AB(i,:)  = local_ferm_AB;

end


%% ================================================================
%  8. Common filter for metabolites and fluxes
%  ================================================================

% A row is retained only if:
%   1. A satisfies the flux constraints
%   2. B satisfies the flux constraints
%   3. Every AB hybrid satisfies the flux constraints

flux_good = ...
    valid_flux_A & ...
    valid_flux_B & ...
    all(valid_flux_AB, 2);

num_final_valid = sum(flux_good);

fprintf( ...
    'Samples passing all flux constraints: %d/%d\n', ...
    num_final_valid, ...
    num_ode_valid ...
);

fprintf( ...
    'Final fraction of original Sobol samples retained: %.2f%%\n', ...
    100 * num_final_valid / num_samples ...
);

if num_final_valid < 2
    error( ...
        ['Fewer than two valid Sobol samples remain after applying ', ...
         'the flux constraints.'] ...
    );
end


%% Apply exactly the same filter to every output

steady_YA  = steady_YA(flux_good,:);
steady_YB  = steady_YB(flux_good,:);
steady_YAB = steady_YAB(flux_good,:,:);

tca_flux_A  = tca_flux_A(flux_good);
tca_flux_B  = tca_flux_B(flux_good);
tca_flux_AB = tca_flux_AB(flux_good,:);

ferm_flux_A  = ferm_flux_A(flux_good);
ferm_flux_B  = ferm_flux_B(flux_good);
ferm_flux_AB = ferm_flux_AB(flux_good,:);


%% ================================================================
%  9. Combine metabolite and flux outputs
%  ================================================================

% Convert fluxes to intended units.
% Multiplication by 60 does not change the Sobol indices.

tca_flux_A  = tca_flux_A  * 60;
tca_flux_B  = tca_flux_B  * 60;
tca_flux_AB = tca_flux_AB * 60;

ferm_flux_A  = ferm_flux_A  * 60;
ferm_flux_B  = ferm_flux_B  * 60;
ferm_flux_AB = ferm_flux_AB * 60;


% Dimensions of Y_A and Y_B:
% sample x output

Y_A = [
    steady_YA, ...
    tca_flux_A, ...
    ferm_flux_A
];

Y_B = [
    steady_YB, ...
    tca_flux_B, ...
    ferm_flux_B
];


% Dimensions of Y_AB:
% sample x parameter x output

Y_AB = nan( ...
    num_final_valid, ...
    num_params, ...
    num_total_outputs ...
);

Y_AB(:,:,1:num_metabolite_outputs) = steady_YAB;

Y_AB(:,:,14) = reshape( ...
    tca_flux_AB, ...
    num_final_valid, ...
    num_params, ...
    1 ...
);

Y_AB(:,:,15) = reshape( ...
    ferm_flux_AB, ...
    num_final_valid, ...
    num_params, ...
    1 ...
);


%% ================================================================
%  10. Calculate Sobol indices for all 15 outputs
%  ================================================================

[S_all, ST_all, output_variance] = calculate_sobol_indices( ...
    Y_A, ...
    Y_B, ...
    Y_AB ...
);

output_names = { ...
    'Pyruvate', ...
    'Oxaloacetate', ...
    'Acetyl-CoA', ...
    'Citrate', ...
    'Cis-aconitate', ...
    'Isocitrate', ...
    '\alpha-Ketoglutarate', ...
    'Succinyl-CoA', ...
    'Succinate', ...
    'Fumarate', ...
    'Acetyl-P', ...
    'NAD', ...
    'NADH', ...
    'TCA flux', ...
    'Fermentation flux' ...
};


%% ================================================================
%  11. Sort parameters using maximum total-order index
%  ================================================================

[mean_total_order, sort_idx] = sort( ...
    median(ST_all, 2), ...
    'descend' ...
);

S_sorted = S_all(sort_idx,:);
ST_sorted = ST_all(sort_idx,:);
fields_sorted = fields(sort_idx);



%% Create paper-ready parameter labels

parameter_label_map = containers.Map( ...
    { ...
        'nad_pool', ...
        'enzyme_sd', ...
        'enzyme_ack', ...
        'alpha', ...
        'enzyme_id', ...
        'enzyme_akd', ...
        'ki_nadh_pdh', ...
        'enzyme_pdh', ...
        'enzyme_pk', ...
        'enzyme_fumarase', ...
        'enzyme_pta', ...
        'k', ...
        'enzyme_aco', ...
        'enzyme_cs', ...
        'ki_nadh_cs', ...
        'enzyme_md', ...
        'enzyme_scs' ...
    }, ...
    { ...
        'NAD(H) pool', ...
        'E_{SDH}', ...
        'E_{ACKA}', ...
        '\alpha', ...
        'E_{IDH}', ...
        'E_{AKGDH}', ...
        'K_{i, NADH}^{PDH}', ...
        'E_{PDH}', ...
        'E_{PK}', ...
        'E_{FUM}', ...
        'E_{PTA}', ...
        'k', ...
        'E_{ACO}', ...
        'E_{CS}', ...
        'K_{i, NADH}^{CS}', ...
        'E_{MDH}', ...
        'E_{SCS}' ...
    } ...
);

parameter_display_names = cell(num_params, 1);

for i = 1:num_params

    field_name = fields_sorted{i};

    if isKey(parameter_label_map, field_name)
        parameter_display_names{i} = ...
            parameter_label_map(field_name);
    else
        % Fallback for any parameter not explicitly listed above
        parameter_display_names{i} = ...
            strrep(field_name, '_', '\_');
    end

end


%% ================================================================
%  12. Colour limits
%  ================================================================

colour_min = 0;
colour_max = 0.545454545455;   % Halve the unlabelled gap above the 0.5 tick.

% Always include zero
colour_min = min(0, colour_min);
colour_max = max(0, colour_max);

if colour_max <= colour_min
    colour_max = colour_min + 1;
end


%% ================================================================
%  13. Plot total-order Sobol heatmap
%  ================================================================

[fig1, ax1] = make_journal_figure(24, 24);

[h2, colourbar2] = plot_sobol_heatmap( ...
    fig1, ax1, ST_sorted, output_names, parameter_display_names, ...
    [colour_min, colour_max] ...
);


%% ================================================================
%  13.2. Plot first-order Sobol heatmap
%  ================================================================

[fig2, ax2] = make_journal_figure(28, 28);

[h3, colourbar3] = plot_sobol_heatmap( ...
    fig2, ax2, S_sorted, output_names, parameter_display_names, ...
    [colour_min, colour_max] ...
);


%% ================================================================
%  14. Create long-format results table
%  ================================================================

%% Create long-format results table

[parameter_index, output_index] = ndgrid( ...
    1:num_params, ...
    1:num_total_outputs ...
);

% Explicitly convert every table variable to a column vector
parameter_column = string(fields(parameter_index(:)));
parameter_column = parameter_column(:);

output_column = string(output_names(output_index(:)));
output_column = output_column(:);

first_order_column = S_all(:);
total_order_column = ST_all(:);

% Check that all columns have identical lengths
assert( ...
    numel(parameter_column) == numel(output_column) && ...
    numel(output_column) == numel(first_order_column) && ...
    numel(first_order_column) == numel(total_order_column), ...
    'The Sobol result columns have inconsistent lengths.' ...
);

sobol_results = table( ...
    parameter_column, ...
    output_column, ...
    first_order_column, ...
    total_order_column, ...
    'VariableNames', { ...
        'Parameter', ...
        'Output', ...
        'FirstOrder', ...
        'TotalOrder' ...
    } ...
);

% Sort from largest to smallest total-order index
sobol_results = sortrows( ...
    sobol_results, ...
    'TotalOrder', ...
    'descend' ...
);

fprintf('\nTop 30 parameter-output combinations:\n\n');

disp(sobol_results(1:min(30, height(sobol_results)), :));

%% ================================================================
%  Local functions
%  All functions must remain at the end of the script
%  ================================================================

function [image_handle, colour_bar] = plot_sobol_heatmap( ...
    fig, ax, indices, output_names, parameter_labels, colour_limits ...
)
% Use regular axes so the colour bar can be positioned independently.

    image_handle = imagesc(ax, indices, colour_limits);

    % Keep the original blues through 0.5, then deepen towards the top.
    deep_blue = [0, 0.4470, 0.7410];
    pale_blue = 0.9 + 0.1 * deep_blue;
    darkest_blue = [0, 0.25, 0.55];
    colour_values = linspace(colour_limits(1), colour_limits(2), 256)';
    blue_map = interp1( ...
        [colour_limits(1), 0.5, colour_limits(2)], ...
        [pale_blue; deep_blue; darkest_blue], colour_values ...
    );
    colormap(ax, blue_map);

    set(ax, ...
        'Units', 'normalized', ...
        'FontName', 'Arial', ...
        'FontSize', 15, ...
        'TickLabelInterpreter', 'tex', ...
        'XTick', 1:numel(output_names), ...
        'XTickLabel', output_names, ...
        'XTickLabelRotation', 45, ...
        'YTick', 1:numel(parameter_labels), ...
        'YTickLabel', parameter_labels, ...
        'YDir', 'reverse', ...
        'TickLength', [0, 0], ...
        'Box', 'on', ...
        'LineWidth', 0.8 ...
    );

    xlabel(ax, 'Model output');
    ylabel(ax, 'Parameter');

    colour_bar = colorbar(ax, 'eastoutside');
    colour_bar.Units = 'normalized';
    colour_bar.FontName = 'Arial';
    colour_bar.FontSize = 13.5;
    colour_bar.TickLabelInterpreter = 'tex';
    colour_bar.TickDirection = 'out';
    % The 0.5 tick sits below the top; the darker extension has no label.
    colour_bar.Ticks = 0:0.05:0.5;
    colour_bar.TickLabels = arrayfun( ...
        @(value) sprintf('%g', value), colour_bar.Ticks, ...
        'UniformOutput', false ...
    );

    % Leave a clear gap between the heatmap and its colour bar.
    ax.Position = [0.22, 0.25, 0.64, 0.70];
    colour_bar.Position = [0.89, 0.25, 0.02, 0.70];

    fig.Color = 'w';
end


function [ok, yend] = run_one(p)

    num_outputs = 13;

    yend = nan(1, num_outputs);
    ok   = false;

    tspan = [0, 60 * 60 * 24];

    options = odeset( ...
        'RelTol', 1e-4, ...
        'AbsTol', 1e-6 ...
    );

    y0 = [
        p.pyruvate_initial;
        p.oxaloacetate_initial;
        p.acetyl_coa_initial;
        p.citrate_initial;
        p.cis_aconitate_initial;
        p.isocitrate_initial;
        p.aketo_initial;
        p.succinyl_coa_initial;
        p.succinate_initial;
        p.fumarate_initial;
        p.acetyl_p_initial;
        p.nad_pool * (1 - p.nadh_ratio);  % NAD
        p.nad_pool * p.nadh_ratio        % NADH
    ];

    try

        [t, y] = ode15s( ...
            @(t,y) function_fsolve(y,p), ...
            tspan, ...
            y0, ...
            options ...
        );

        if isempty(t) || isempty(y)
            return;
        end

        reached_final_time = ...
            t(end) >= 0.9 * tspan(end);

        finite_solution = ...
            all(isfinite(y(end,:)));

        bounded_solution = ...
            all(abs(y(end,:)) < 1e12);

        correct_output_size = ...
            size(y,2) == num_outputs;

        ok = ...
            numel(t) > 1 && ...
            reached_final_time && ...
            finite_solution && ...
            bounded_solution && ...
            correct_output_size;

        if ok
            yend = reshape(y(end,:), 1, num_outputs);
        end

    catch
        ok   = false;
        yend = nan(1, num_outputs);
    end

end


function [is_valid, tca_flux, ferm_flux] = ...
    extract_flux_outputs(y, p)

    tca_flux  = nan;
    ferm_flux = nan;
    is_valid  = false;

    % Flux consistency criterion:
    % allow 1% relative difference, with a small absolute tolerance
    absolute_flux_tolerance = 1e-6;
    relative_flux_tolerance = 1e-2;   % 1%

    try

        fluxes = compute_flux(y, p);

        if isempty(fluxes)
            return;
        end

        % Use the final row if compute_flux returns multiple time points
        fluxes = fluxes(end,:);

        % Force into a row vector
        fluxes = fluxes(:).';

        if numel(fluxes) < 12 || any(~isfinite(fluxes(1:12)))
            return;
        end

        tca_reaction_fluxes = fluxes(3:11);

        %% Check whether TCA reaction fluxes are consistent

        flux_range = ...
            max(tca_reaction_fluxes) - ...
            min(tca_reaction_fluxes);

        flux_scale = max(abs(tca_reaction_fluxes));

        allowed_difference = ...
            absolute_flux_tolerance + ...
            relative_flux_tolerance * flux_scale;

        tca_consistent = ...
            flux_range <= allowed_difference;


        %% Allow zero flux, but reject clearly reverse TCA flux

        not_substantially_reverse = ...
            median(tca_reaction_fluxes) >= ...
            -absolute_flux_tolerance;


        %% Final validity criterion

        is_valid = ...
            tca_consistent && ...
            not_substantially_reverse;


        %% Extract outputs

        if is_valid

            % Use the median of the consistent TCA reaction fluxes
            tca_flux = median(tca_reaction_fluxes);

            ferm_flux = fluxes(12);

        end

    catch

        is_valid  = false;
        tca_flux  = nan;
        ferm_flux = nan;

    end

end


function [S, ST, D] = calculate_sobol_indices( ...
    Y_A, ...
    Y_B, ...
    Y_AB ...
)

    num_samples = size(Y_A, 1);
    num_outputs = size(Y_A, 2);
    num_params  = size(Y_AB, 2);

    assert( ...
        isequal(size(Y_A), size(Y_B)), ...
        'Y_A and Y_B must have identical dimensions.' ...
    );

    assert( ...
        size(Y_AB,1) == num_samples, ...
        'Y_AB has the wrong number of samples.' ...
    );

    assert( ...
        size(Y_AB,3) == num_outputs, ...
        'Y_AB has the wrong number of outputs.' ...
    );

    if any(~isfinite(Y_A), 'all') || ...
       any(~isfinite(Y_B), 'all') || ...
       any(~isfinite(Y_AB), 'all')

        error( ...
            ['Non-finite values remain after filtering. ', ...
             'Check the common validity mask.'] ...
        );
    end

    % Variance of each output using combined A and B evaluations
    D = var([Y_A; Y_B], 0, 1);

    if any(~isfinite(D)) || any(D <= 0)

        bad_outputs = find(~isfinite(D) | D <= 0);

        error( ...
            ['The following output columns have zero or invalid ', ...
             'variance: %s'], ...
            mat2str(bad_outputs) ...
        );

    end

    S  = nan(num_params, num_outputs);
    ST = nan(num_params, num_outputs);

    for h = 1:num_params

        Y_AB_h = reshape( ...
            Y_AB(:,h,:), ...
            num_samples, ...
            num_outputs ...
        );

        % First-order estimator
        S(h,:) = mean( ...
            Y_B .* (Y_AB_h - Y_A), ...
            1 ...
        ) ./ D;

        % Total-order estimator
        ST(h,:) = mean( ...
            (Y_A - Y_AB_h).^2, ...
            1 ...
        ) ./ (2 .* D);

    end

end
