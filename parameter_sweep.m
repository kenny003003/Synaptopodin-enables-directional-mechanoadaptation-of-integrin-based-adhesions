clear; close all; clc;

%% === PARAMETER SWEEP FOR DETACHMENT RATE ===
% Sweeps through F_inertial values, simulates 1000 cells per condition
% Compares two c values and mixed aspect ratio populations

fprintf('=== Detachment Rate Parameter Sweep ===\n');
fprintf('Starting simulation...\n\n');

%% --- Fixed Parameters ---
r0 = 0.5;        % Actin cap radius
h = 0.01;         % Height above substrate
k0 = 10000;          % Stiffness of fiber with length h
F0 = 0.0;        % Actomyosin contractile force

% Off-rate parameters
k_off0 = 1;
N0 = 1;
K_s = 1;
K_0 = 1;
n_exp = 1;
m_exp = 2;

% Simulation parameters
N = 100;         % Number of fibers per cell
dt = 0.1;      % Time step
n_timesteps = 25;  % Simulation duration
max_iter = 50;   
tol = 1e-6;      

% Fiber distribution
POLARIZATION = 0;
USE_RANDOM = true;
JITTER_AMOUNT = 0.05;

%% --- Sweep Parameters ---
F_inertial_values = 15 * [0.067, 0.23, 1];  % [0.67, 2.3, 10]
c_values = [1];  % Two off-rate parameters to compare
n_cells = 100;  % Number of cells to simulate per condition

% Aspect ratios
AR_low = 1;   % Low aspect ratio (nearly circular)
AR_high = 8;    % High aspect ratio (elongated)
a_fixed = 1;

% Population composition
fraction_AR_low = 1;  % Fraction of cells with low AR (0 to 1)
                        % 0.5 = 50% low AR, 50% high AR
                        % 0.0 = all cells are high AR
                        % 1.0 = all cells are low AR

% Detachment criterion
detachment_threshold = 0.1;  % Cell is detached if < 5% fibers remain

fprintf('Simulation parameters:\n');
fprintf('  Forces: [%.3f, %.3f, %.3f]\n', F_inertial_values);
fprintf('  c values: [%d, %d]\n', c_values);
fprintf('  Cells per condition: %d\n', n_cells);
fprintf('  AR low: %.1f, AR high: %.1f\n', AR_low, AR_high);
fprintf('  Population: %.0f%% low AR, %.0f%% high AR\n', ...
        fraction_AR_low*100, (1-fraction_AR_low)*100);
fprintf('  Detachment threshold: < %.0f%% fibers (<%d fibers for N=%d)\n', ...
        detachment_threshold*100, round(detachment_threshold*N), N);
fprintf('  Simulation time: %.2f seconds\n', n_timesteps*dt);
fprintf('  Polarization: %.2f\n\n', POLARIZATION);

%% --- Initialize storage ---
n_forces = length(F_inertial_values);
n_c_values = length(c_values);

% Store results: detachment_rate(c_index, force_index)
% Detachment rate = fraction of cells with < 5% fibers remaining
detachment_rate = zeros(n_c_values, n_forces);
n_detached_cells = zeros(n_c_values, n_forces);  % Count of detached cells

% Determine AR for each cell (fixed across all conditions)
n_AR_low = round(n_cells * fraction_AR_low);
n_AR_high = n_cells - n_AR_low;
cell_AR = [AR_low * ones(1, n_AR_low), AR_high * ones(1, n_AR_high)];
cell_AR = cell_AR(randperm(n_cells));  % Shuffle

%% --- Run parameter sweep ---
for c_idx = 1:n_c_values
    c = c_values(c_idx);
    fprintf('Testing c = %d:\n', c);
    
    for f_idx = 1:n_forces
        F_inertial = F_inertial_values(f_idx);
        fprintf('  F_inertial = %.3f: ', F_inertial);
        
        % Store detachment status for all cells (binary: detached or not)
        cell_is_detached = false(1, n_cells);
        
        % Simulate each cell
        for cell_idx = 1:n_cells
            % Get aspect ratio for this cell
            AR = cell_AR(cell_idx);
            b = AR * a_fixed;
            
            % Simulate this cell
            [bonded_hist, ~] = simulate_cell(a_fixed, b, r0, h, k0, F0, c, ...
                F_inertial, N, dt, n_timesteps, max_iter, tol, ...
                k_off0, N0, K_s, K_0, n_exp, m_exp, ...
                POLARIZATION, USE_RANDOM, [], JITTER_AMOUNT);
            
            % Check if cell is detached (< 5% fibers remain)
            final_bonded_fraction = sum(bonded_hist(end,:)) / N;
            if final_bonded_fraction < detachment_threshold
                cell_is_detached(cell_idx) = true;
            end
            
            % Progress indicator
            if mod(cell_idx, 100) == 0
                fprintf('.');
            end
        end
        
        % Calculate detachment rate (fraction of cells that detached)
        n_detached = sum(cell_is_detached);
        detach_rate = n_detached / n_cells;
        
        detachment_rate(c_idx, f_idx) = detach_rate;
        n_detached_cells(c_idx, f_idx) = n_detached;
        
        fprintf(' Done! Detachment rate: %.3f (%d/%d cells)\n', ...
                detach_rate, n_detached, n_cells);
    end
    fprintf('\n');
end

fprintf('Parameter sweep complete!\n\n');

%% --- Plot results ---
figure('Position', [100 100 900 600]);
hold on;

% Color scheme
colors = [0.2 0.4 0.8; 0.8 0.2 0.2];  % Blue, Red
markers = {'o-', 's-'};
line_width = 2.5;
marker_size = 10;

% Plot each c value
for c_idx = 1:n_c_values
    plot(F_inertial_values, detachment_rate(c_idx, :), markers{c_idx}, ...
         'Color', colors(c_idx, :), ...
         'LineWidth', line_width, ...
         'MarkerSize', marker_size, ...
         'MarkerFaceColor', colors(c_idx, :), ...
         'DisplayName', sprintf('c = %d', c_values(c_idx)));
end

% Formatting
xlabel('Applied Force (F_{inertial})', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Detachment Rate (fraction of cells)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Detachment Rate vs Applied Force\n(%d cells, %.0f%% AR=%.1f, %.0f%% AR=%.1f, Threshold=%.0f%%)', ...
      n_cells, fraction_AR_low*100, AR_low, (1-fraction_AR_low)*100, AR_high, detachment_threshold*100), ...
      'FontSize', 15, 'FontWeight', 'bold');

legend('Location', 'northwest', 'FontSize', 12);
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

% Set x-axis to log scale if span is large
if max(F_inertial_values) / min(F_inertial_values) > 10
    set(gca, 'XScale', 'log');
end

ylim([0, 1.1]);
hold off;

%% --- Print summary table ---
fprintf('=== Summary Table ===\n');
fprintf('Cell is detached if < %.0f%% fibers remain (< %d fibers for N=%d)\n\n', ...
        detachment_threshold*100, round(detachment_threshold*N), N);
fprintf('Force    ');
for c_idx = 1:n_c_values
    fprintf('| c=%4d         ', c_values(c_idx));
end
fprintf('\n');
fprintf('---------');
for c_idx = 1:n_c_values
    fprintf('+---------------');
end
fprintf('\n');

for f_idx = 1:n_forces
    fprintf('%.3f    ', F_inertial_values(f_idx));
    for c_idx = 1:n_c_values
        fprintf('| %.3f (%4d/%4d)', ...
                detachment_rate(c_idx, f_idx), ...
                n_detached_cells(c_idx, f_idx), ...
                n_cells);
    end
    fprintf('\n');
end
fprintf('\n');

%% --- Function to simulate a single cell ---
function [bonded_history, delta_history] = simulate_cell(a, b, r0, h, k0, F0, c, ...
    F_inertial, N, dt, n_timesteps, max_iter, tol, ...
    k_off0, N0, K_s, K_0, n_exp, m_exp, ...
    polarization, use_random, random_seed, jitter_amount)
    
    % Setup fibers with polarized distribution
    theta_i = generate_fiber_angles(N, polarization, use_random, random_seed, jitter_amount);
    r_theta = a ./ sqrt(cos(theta_i).^2 + (a/b)^2 * sin(theta_i).^2);
    L0 = sqrt((r_theta - r0).^2 + h^2);
    k_i = k0 * h ./ L0;
    
    % Initialize storage
    bonded_history = true(n_timesteps, N);
    delta_history = zeros(n_timesteps, 1);
    
    % Initialize bonded status
    bonded = true(1, N);
    delta_prev = 0;
    
    % Time evolution
    for j = 1:n_timesteps
        if sum(bonded) < 5
            % Cell fully detached
            bonded_history(j:end, :) = false;
            delta_history(j:end) = NaN;
            break;
        end
        
        % Solve for displacement
        delta_j = delta_prev;
        for iter = 1:max_iter
            phi_i = atan((h + delta_j) ./ (r_theta - r0));
            L_f = sqrt((r_theta - r0).^2 + (h + delta_j)^2);
            F_i = F0 + k_i .* h .* (L_f ./ L0 - 1);
            F_sum = sum(F_i(bonded) .* sin(phi_i(bonded)));
            residual = F_sum - F_inertial;
            
            if abs(residual) < tol
                break;
            end
            
            step_size = 0.1;
            delta_j = delta_j - step_size * residual;
            if delta_j < 0
                delta_j = 0;
                break;
            end
        end
        
        delta_prev = delta_j;
        
        % Calculate final forces for this timestep
        phi_i = atan((h + delta_j) ./ (r_theta - r0));
        L_f = sqrt((r_theta - r0).^2 + (h + delta_j)^2);
        F_i = F0 + k_i .* h .* (L_f ./ L0 - 1);
        N_i = F_i .* sin(phi_i);
        S_i = F_i .* cos(phi_i);
        
        % Off-rate calculation
        k_off = zeros(1, N);
        k_off(bonded) = k_off0 * exp(N_i(bonded)/N0 .* ...
            (K_s^n_exp ./ (K_s^n_exp + S_i(bonded).^n_exp)) .* ...
            (1 + c * (S_i(bonded).^m_exp ./ (K_0^m_exp + S_i(bonded).^m_exp))));
        
        % Debonding probability
        % prob_debond = k_off * dt;
        prob_debond = 1-exp(-k_off*dt);
        
        % Store this timestep
        bonded_history(j, :) = bonded;
        delta_history(j) = delta_j;
        
        % Debonding
        random_chance = rand(1, N);
        newly_broken = bonded & (prob_debond >= random_chance);
        bonded(newly_broken) = false;
    end
end

%% --- Function to generate polarized fiber distribution ---
function theta_i = generate_fiber_angles(N, polarization, use_random, random_seed, jitter_amount)
    % Generate fiber angles with optional polarization toward long axis
    
    % Set random seed if specified
    if ~isempty(random_seed)
        rng(random_seed);
    end
    
    if polarization < 0.05
        % Uniform distribution
        if use_random
            theta_i = sort(rand(1, N) * 2 * pi);
        else
            theta_i = linspace(0, 2*pi*(N-1)/N, N);
        end
    else
        % Generate non-uniform distribution concentrated at long axis
        power = 1 + polarization * 5;
        
        theta_fine = linspace(0, 2*pi, 10000);
        density = (sin(theta_fine).^2).^power;
        density = density + 1e-10;
        density = density / sum(density);
        cdf = cumsum(density);
        cdf = cdf / cdf(end);
        
        [cdf_unique, unique_idx] = unique(cdf);
        theta_unique = theta_fine(unique_idx);
        
        if use_random
            u = rand(1, N);
            u = sort(u);
        else
            u = linspace(0, 1, N);
        end
        
        u = max(min(u, cdf_unique(end)), cdf_unique(1));
        theta_i = interp1(cdf_unique, theta_unique, u, 'linear', 'extrap');
        
        if use_random || jitter_amount > 0
            jitter = jitter_amount * randn(size(theta_i));
            theta_i = theta_i + jitter;
        end
        
        theta_i = mod(theta_i, 2*pi);
    end
end