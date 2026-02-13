clear; close all; clc;

%% === PARAMETERS ===
% Sweep Settings
F_BASE = 4;
FORCE_RATIOS = [1, 3.5, 15];
C_VALUES = [2.5];%%%0.7 for control%%%%1.3 for KO
AR_VALUES = [2];

% Simulation Control
N_TRIALS = 1000;            % Cells per condition
N_FIBERS = 100;             % Fibers per cell
DETACH_THRESHOLD = 0.05;    % <5% fibers = detached
MAX_TIMESTEPS = 250;
DT = 0.01;

% Physical Constants
r0 = 0.3; h = 0.05; k0 = 10000000; F0 = 0.0;
k_off0 = 1; N0 = 1; K_s = 1; K_0 = 1; n_exp = 1; m_exp = 2;

%% === SIMULATION LOOP ===
% Storage: [C_idx, AR_idx, Force_idx]
results_detach = zeros(length(C_VALUES), length(AR_VALUES), length(FORCE_RATIOS));
results_conv   = zeros(length(C_VALUES), length(AR_VALUES), length(FORCE_RATIOS));

% Check for parallel pool
if isempty(gcp('nocreate'))
    parpool; % Start parallel pool if not already running
end

fprintf('Starting simulation (%d trials per condition)...\n', N_TRIALS);
t_start = tic;

for c_idx = 1:length(C_VALUES)
    C = C_VALUES(c_idx);
    
    for ar_idx = 1:length(AR_VALUES)
        AR = AR_VALUES(ar_idx);
        
        % Set Polarization: 0.7 for elongated (AR>1.5), 0.0 for circular
        if AR > 1.5
            polarization = 0.7;
        else
            polarization = 0.0;
        end
        
        % --- MODIFIED GEOMETRY (Constant Area) ---
        radius_base = 1.0;
        a_fixed = radius_base / sqrt(AR);
        b_axis = a_fixed * AR;
        % -----------------------------------------
        
        for f_idx = 1:length(FORCE_RATIOS)
            F_ratio = FORCE_RATIOS(f_idx);
            F_inertial = F_ratio * F_BASE;
            
            % Initialize reduction variables for parfor
            detached_count = 0;
            converged_trials_count = 0;
            
            % --- PARALLEL TRIALS ---
            parfor trial = 1:N_TRIALS
                % 1. Generate Fibers
                % (RNG is handled automatically by parfor workers)
                theta_i = generate_fiber_angles(N_FIBERS, polarization);
                r_theta = a_fixed ./ sqrt(cos(theta_i).^2 + (a_fixed/b_axis)^2 * sin(theta_i).^2);
                L0 = sqrt((r_theta - r0).^2 + h^2);
                k_i = k0 * h ./ L0;
  
                % 2. Initialize Local Variables
                bonded = true(1, N_FIBERS);
                delta = 0;
                is_detached = false;
                all_converged = true;
                
                % 3. Time Steps
                for t = 1:MAX_TIMESTEPS
                    if sum(bonded) < (N_FIBERS * DETACH_THRESHOLD)
                        is_detached = true;
                        break;
                    end
                    
                    % Equilibrium Solver
                    [delta, converged] = solve_equilibrium(delta, r_theta, r0, h, k_i, L0, F0, bonded, F_inertial);
                    
                    if ~converged
                        all_converged = false;
                    end
                    
                    % Forces & Kinetics
                    phi_i = atan((h + delta) ./ (r_theta - r0));
                    L_f = sqrt((r_theta - r0).^2 + (h + delta).^2);
                    F_i = F0 + k_i .* h .* (L_f ./ L0 - 1);
                    
                    N_i = F_i .* sin(phi_i);
                    S_i = F_i .* cos(phi_i);
                    
                    % Bond Logic
                    k_off = zeros(1, N_FIBERS);
                    mask = bonded;
                    if any(mask)
                        term1 = (K_s^n_exp) ./ (K_s^n_exp + S_i(mask).^n_exp);
                        term2 = 1 + C .* (S_i(mask).^m_exp ./ (K_0^m_exp + S_i(mask).^m_exp));
                        k_off(mask) = k_off0 .* exp(N_i(mask)./N0 .* term1 .* term2);
                    end
                    
                    prob = 1 - exp(-k_off * DT);
                    bonded(bonded & (rand(1, N_FIBERS) <= prob)) = false;
                end
                
                % Accumulate Results (Reduction)
                if is_detached
                    detached_count = detached_count + 1;
                end
                if all_converged
                    converged_trials_count = converged_trials_count + 1;
                end
            end
            % --- END PARFOR ---
            
            % Store Results
            results_detach(c_idx, ar_idx, f_idx) = detached_count / N_TRIALS;
            results_conv(c_idx, ar_idx, f_idx) = converged_trials_count / N_TRIALS;
            
            fprintf('C=%.1f, AR=%.1f, F=%.1fx -> Detach: %.2f, Conv: %.2f\n', ...
                    C, AR, F_ratio, results_detach(c_idx, ar_idx, f_idx), results_conv(c_idx, ar_idx, f_idx));
        end
    end
end

toc(t_start);

%% === PLOTTING ===
figure('Position', [100 100 1200 500]);

% Styles matching the Python plot
colors = {'b', 'r'};       % Blue for C=0.5, Red for C=1.0
markers = {'o', 's'};      % Circle for AR=1, Square for AR=2.2
lines = {'--', '-'};       % Dashed for AR=1, Solid for AR=2.2

% --- Subplot 1: Detachment ---
subplot(1, 2, 1); hold on;
for c_idx = 1:length(C_VALUES)
    for ar_idx = 1:length(AR_VALUES)
        y_data = squeeze(results_detach(c_idx, ar_idx, :));
        label_str = sprintf('C=%.1f, AR=%.1f', C_VALUES(c_idx), AR_VALUES(ar_idx));
        
        plot(FORCE_RATIOS, y_data, ...
             'Color', colors{c_idx}, ...
             'Marker', markers{ar_idx}, ...
             'LineStyle', lines{ar_idx}, ...
             'LineWidth', 1.5, ...
             'DisplayName', label_str);
    end
end
title('Cell Detachment vs Force (Constant Area)');
xlabel('Force Ratio (x 1000 pN)');
ylabel('Detachment Ratio');
ylim([-0.05 1.05]);
grid on;
legend('Location', 'best');

% --- Subplot 2: Convergence ---
subplot(1, 2, 2); hold on;
for c_idx = 1:length(C_VALUES)
    for ar_idx = 1:length(AR_VALUES)
        y_data = squeeze(results_conv(c_idx, ar_idx, :));
        label_str = sprintf('C=%.1f, AR=%.1f', C_VALUES(c_idx), AR_VALUES(ar_idx));
        
        plot(FORCE_RATIOS, y_data, ...
             'Color', colors{c_idx}, ...
             'Marker', markers{ar_idx}, ...
             'LineStyle', lines{ar_idx}, ...
             'LineWidth', 1.5, ...
             'DisplayName', label_str);
    end
end
title('Numerical Convergence vs Force');
xlabel('Force Ratio (x 1000 pN)');
ylabel('Convergence Ratio');
ylim([-0.05 1.05]);
grid on;
legend('Location', 'best');


%% === HELPER FUNCTIONS ===

function theta_i = generate_fiber_angles(N, polarization)
    % Generates fiber angles. 
    % polarization=0 -> Uniform
    % polarization>0 -> Weighted towards pi/2 and 3pi/2
    
    if polarization < 0.05
        theta_i = sort(rand(1, N) * 2 * pi);
    else
        power = 1 + polarization * 5;
        theta_fine = linspace(0, 2*pi, 1000);
        density = (sin(theta_fine).^2).^power + 1e-10;
        density = density / sum(density);
        cdf = cumsum(density);
        cdf = cdf / cdf(end);
        
        % Inverse transform
        [cdf_u, idx] = unique(cdf);
        theta_unique = theta_fine(idx);
        
        u = sort(rand(1, N));
        u = max(min(u, cdf_u(end)), cdf_u(1)); % Clip
        theta_i = interp1(cdf_u, theta_unique, u, 'linear');
        
        % Jitter
        theta_i = mod(theta_i + 0.005 * randn(1, N), 2*pi);
    end
end

function [delta, converged] = solve_equilibrium(delta_init, r_theta, r0, h, k_i, L0, F0, bonded, F_inertial)
    % Newton-Raphson solver with adaptive damping
    delta = delta_init;
    damping = 0.3;
    max_iter = 100;
    tol = 1e-6;
    converged = false;
    
    % Filter bonded only
    r_b = r_theta(bonded);
    k_b = k_i(bonded);
    L0_b = L0(bonded);
    
    if isempty(r_b), converged=true; return; end
    
    for i = 1:max_iter
        phi = atan((h + delta) ./ (r_b - r0));
        L_f = sqrt((r_b - r0).^2 + (h + delta).^2);
        F_i = F0 + k_b .* h .* (L_f ./ L0_b - 1);
        
        residual = sum(F_i .* sin(phi)) - F_inertial;
        
        if abs(residual) < tol
            converged = true; return;
        end
        
        % Derivatives
        dF = k_b .* h ./ L0_b .* ((h + delta) ./ L_f);
        dPhi = (r_b - r0) ./ ((r_b - r0).^2 + (h + delta).^2);
        dR = sum(dF .* sin(phi) + F_i .* cos(phi) .* dPhi);
        
        if abs(dR) < 1e-12, dR = sign(dR)*1e-12; end
        if dR == 0, dR = 1e-12; end
        
        step = -residual / dR;
        
        % Adaptive Damping (Simplified check)
        delta_new = max(0, delta + damping * step);
        
        % Check if residual improves
        phi_new = atan((h + delta_new) ./ (r_b - r0));
        L_new = sqrt((r_b - r0).^2 + (h + delta_new).^2);
        F_new = F0 + k_b .* h .* (L_new ./ L0_b - 1);
        res_new = abs(sum(F_new .* sin(phi_new)) - F_inertial);
        
        if res_new < abs(residual)
            delta = delta_new;
            damping = min(damping * 1.2, 1.0);
        else
            damping = max(damping * 0.5, 0.001);
            % Take conservative step
            delta = max(0, delta + damping * step);
        end
    end
end