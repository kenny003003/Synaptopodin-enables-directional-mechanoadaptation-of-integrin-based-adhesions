clear; close all; clc;

%% === VISUALIZATION MODE SELECTION ===
% Choose what you want to see:
VISUALIZATION_MODE = 'snapshots';  % Options:
                                   % 'snapshots' - Show 4 key timesteps side-by-side
                                   % 'animation' - Create animated movie
                                   % 'single'    - Single timestep (original)

% Animation settings (if VISUALIZATION_MODE = 'animation')
SAVE_ANIMATION = false;  % Set to true to save as .gif or .avi
ANIMATION_FILENAME = 'fiber_evolution.gif';
FRAME_DELAY = 0.1;  % seconds between frames

% Snapshot settings (if VISUALIZATION_MODE = 'snapshots')
SNAPSHOT_TIMES = [1, 50, 100, 150, 200, 250];  % Timesteps to show

%% --- Parameters (Matched to parameter_sweep_v2.m) ---
% Force settings
F_BASE = 5;
FORCE_RATIO = 3.5;           % Example case: Ratio 3.5
F_inertial = FORCE_RATIO * F_BASE; 

% Physical Geometry & Stiffness
r0 = 0.3;        % Actin cap radius (updated from 0.5)
h = 0.05;        % Height above substrate
k0 = 10000000;   % Stiffness (⚠️ VERY HIGH - matched to sweep code)
F0 = 0.0;        % Actomyosin contractile force

% Off-rate parameters (Updated to match sweep)
k_off0 = 1;      % Base off-rate (was 0.01)
N0 = 1;
K_s = 1;
K_0 = 1;
n_exp = 1;
m_exp = 2;
c = 0.5;         % Off-rate parameter (C_VALUES from sweep)

% Simulation parameters
N = 100;         % Number of fibers
dt = 0.01;       % Time step (matched to sweep)
n_timesteps = 250; % Duration (matched to sweep)
max_iter = 100;  % Solver iterations
tol = 1e-6;      % Solver tolerance

% ========== CONVERGENCE SETTINGS ==========
USE_NEWTON_RAPHSON = true;   % Essential for k0=1e7
USE_ADAPTIVE_DAMPING = true; % Essential for stability
INITIAL_DAMPING = 0.3;
MIN_DAMPING = 0.001;
MAX_DAMPING = 1.0;
AUTO_ADJUST_K0 = true;       % Allow auto-reduction if physics breaks
K0_REDUCTION_FACTOR = 0.5;   
VERBOSE_CONVERGENCE = false; 

% ========== FIBER POLARIZATION CONTROL ==========
% Controls how fibers are distributed around ellipse
POLARIZATION1 = 0;    % Circular case (Uniform)
POLARIZATION2 = 1;  % Elongated case (Polarized)

% ========== RANDOMNESS CONTROL ==========
USE_RANDOM = true;
RANDOM_SEED = [];    % Empty = random every time
JITTER_AMOUNT = 0.005; 

% ========== CONVERGENCE MONITORING ==========
CHECK_CONVERGENCE = true;  % Monitor and report convergence issues
PLOT_CONVERGENCE = true;   % Create convergence diagnostic plots

% Two aspect ratios to compare
AR1 = 1.0;    % Circular
AR2 = 5;    % Elongated

fprintf('=== Fiber Polarization & Convergence Analysis ===\n');
fprintf('Mode: %s\n', VISUALIZATION_MODE);
fprintf('Force: %.1f pN (Ratio %.1f * Base %.0f)\n', F_inertial, FORCE_RATIO, F_BASE);
fprintf('Stiffness k0: %.1e (High Stiffness)\n', k0);

%% --- Simulate both aspect ratios ---

% Solver Configuration
solver_params = struct('use_newton', USE_NEWTON_RAPHSON, ...
                       'use_adaptive', USE_ADAPTIVE_DAMPING, ...
                       'init_damp', INITIAL_DAMPING, ...
                       'min_damp', MIN_DAMPING, ...
                       'max_damp', MAX_DAMPING, ...
                       'damp_inc', 1.2, ...
                       'damp_dec', 0.5, ...
                       'auto_adjust_k0', AUTO_ADJUST_K0, ...
                       'k0_reduction', K0_REDUCTION_FACTOR, ...
                       'verbose', VERBOSE_CONVERGENCE);

% --- SIMULATION 1: AR = 1.0 (Circular) ---
fprintf('Simulating AR %.1f (Constant Area)...\n', AR1);
% Constant Area Geometry Calculation
radius_base = 1.0;
a1 = radius_base / sqrt(AR1);
b1 = a1 * AR1;

[bonded_hist1, F_hist1, N_hist1, S_hist1, theta_i1, r_theta1, delta_hist1, prob_debond_hist1, conv_stats1, k0_used1] = ...
    simulate_full_history(a1, b1, r0, h, k0, F0, c, F_inertial, N, dt, ...
                         n_timesteps, max_iter, tol, k_off0, N0, K_s, K_0, n_exp, m_exp, ...
                         POLARIZATION1, USE_RANDOM, RANDOM_SEED, JITTER_AMOUNT, CHECK_CONVERGENCE, solver_params);

% --- SIMULATION 2: AR = 2.2 (Elongated) ---
fprintf('Simulating AR %.1f (Constant Area)...\n', AR2);
% Constant Area Geometry Calculation
a2 = radius_base / sqrt(AR2);
b2 = a2 * AR2;

[bonded_hist2, F_hist2, N_hist2, S_hist2, theta_i2, r_theta2, delta_hist2, prob_debond_hist2, conv_stats2, k0_used2] = ...
    simulate_full_history(a2, b2, r0, h, k0, F0, c, F_inertial, N, dt, ...
                         n_timesteps, max_iter, tol, k_off0, N0, K_s, K_0, n_exp, m_exp, ...
                         POLARIZATION2, USE_RANDOM, RANDOM_SEED, JITTER_AMOUNT, CHECK_CONVERGENCE, solver_params);

fprintf('Simulation complete!\n');

%% --- Report convergence statistics ---
if CHECK_CONVERGENCE
    report_convergence(AR1, conv_stats1, k0, k0_used1, max_iter, tol);
    report_convergence(AR2, conv_stats2, k0, k0_used2, max_iter, tol);
end

fprintf('\nCreating visualization...\n\n');

%% --- Execute selected visualization mode ---
switch VISUALIZATION_MODE
    case 'single'
        %% SINGLE TIMESTEP MODE
        timestep = 1;
        figure('Position', [50 50 1600 700]);
        
        ax1 = subplot(1, 2, 1);
        plot_timestep(ax1, a1, b1, r0, theta_i1, r_theta1, bonded_hist1(timestep,:), prob_debond_hist1(timestep,:), AR1, timestep, dt);
        colorbar(ax1);
        
        ax2 = subplot(1, 2, 2);
        plot_timestep(ax2, a2, b2, r0, theta_i2, r_theta2, bonded_hist2(timestep,:), prob_debond_hist2(timestep,:), AR2, timestep, dt);
        colorbar(ax2);
        
        sgtitle(sprintf('Fiber State (F=%.1f, k_0=%.1e)', F_inertial, k0), 'FontSize', 16);
        
    case 'snapshots'
        %% SNAPSHOT MODE
        n_snapshots = length(SNAPSHOT_TIMES);
        figure('Position', [50 50 1800 900]);
        
        for i = 1:n_snapshots
            t_idx = SNAPSHOT_TIMES(i);
            if t_idx > n_timesteps, continue; end
            
            % AR 1
            ax = subplot(2, n_snapshots, i);
            plot_timestep(ax, a1, b1, r0, theta_i1, r_theta1, bonded_hist1(t_idx,:), prob_debond_hist1(t_idx,:), AR1, t_idx, dt);
            if i == n_snapshots, colorbar(ax); end
            
            % AR 2
            ax = subplot(2, n_snapshots, i + n_snapshots);
            plot_timestep(ax, a2, b2, r0, theta_i2, r_theta2, bonded_hist2(t_idx,:), prob_debond_hist2(t_idx,:), AR2, t_idx, dt);
            if i == n_snapshots, colorbar(ax); end
        end
        
        sgtitle(sprintf('Evolution Snapshots (F=%.1f, k_0=%.1e)', F_inertial, k0), 'FontSize', 16);
        
    case 'animation'
        %% ANIMATION MODE
        fig = figure('Position', [50 50 1600 700]);
        % (Animation logic omitted for brevity - same as original but with new geometry variables a1, b1, a2, b2)
        for timestep = 1:5:n_timesteps % Skip frames for speed
            clf;
            ax1 = subplot(1, 2, 1);
            plot_timestep(ax1, a1, b1, r0, theta_i1, r_theta1, bonded_hist1(timestep,:), prob_debond_hist1(timestep,:), AR1, timestep, dt);
            
            ax2 = subplot(1, 2, 2);
            plot_timestep(ax2, a2, b2, r0, theta_i2, r_theta2, bonded_hist2(timestep,:), prob_debond_hist2(timestep,:), AR2, timestep, dt);
            
            sgtitle(sprintf('Time: %.2fs', timestep*dt));
            drawnow;
        end
end

%% --- Convergence Diagnostic Plots ---
if CHECK_CONVERGENCE && PLOT_CONVERGENCE
    plot_convergence_diagnostics(conv_stats1, conv_stats2, AR1, AR2, max_iter, tol, k0, F_inertial, n_timesteps, bonded_hist1, bonded_hist2, delta_hist1, delta_hist2, N);
end

%% ========================================================================
%%                      HELPER FUNCTIONS
%% ========================================================================

function report_convergence(AR, stats, k0_orig, k0_used, max_iter, tol)
    fprintf('\nAR %.1f Convergence:\n', AR);
    if k0_used ~= k0_orig
        fprintf('  ⚠️  k0 auto-adjusted: %.1e → %.1e\n', k0_orig, k0_used);
    end
    fprintf('  Avg Iterations: %.1f / %d\n', stats.avg_iter, max_iter);
    fprintf('  Failure Rate: %.1f%%\n', stats.fail_rate * 100);
    fprintf('  Avg Residual: %.2e (tol=%.2e)\n', stats.avg_residual, tol);
end

function plot_convergence_diagnostics(cs1, cs2, AR1, AR2, max_iter, tol, k0, F, n_steps, bh1, bh2, dh1, dh2, N)
    figure('Position', [100 50 1400 900]);
    
    % 1. Iteration History
    subplot(2, 2, 1);
    plot(cs1.iter_history, 'b-', 'DisplayName', sprintf('AR %.1f', AR1)); hold on;
    plot(cs2.iter_history, 'r-', 'DisplayName', sprintf('AR %.1f', AR2));
    yline(max_iter, 'k--', 'Limit');
    title('Solver Iterations per Timestep'); ylabel('Iterations'); legend; grid on;
    
    % 2. Residual History
    subplot(2, 2, 2);
    semilogy(cs1.residual_history, 'b-', 'DisplayName', sprintf('AR %.1f', AR1)); hold on;
    semilogy(cs2.residual_history, 'r-', 'DisplayName', sprintf('AR %.1f', AR2));
    yline(tol, 'k--', 'Tolerance');
    title('Final Residual (Log Scale)'); ylabel('Residual'); grid on;
    
    % 3. Bonded Count
    subplot(2, 2, 3);
    plot(sum(bh1, 2), 'b-', 'LineWidth', 1.5); hold on;
    plot(sum(bh2, 2), 'r-', 'LineWidth', 1.5);
    title('Active Bonds'); ylabel('Count'); ylim([0 N]); grid on;
    
    % 4. Displacement
    subplot(2, 2, 4);
    plot(dh1, 'b-', 'LineWidth', 1.5); hold on;
    plot(dh2, 'r-', 'LineWidth', 1.5);
    title('Cell Displacement'); ylabel('Delta'); grid on;
    
    sgtitle(sprintf('Convergence Diagnostics (k0=%.1e, F=%.1f)', k0, F));
end

function theta_i = generate_fiber_angles(N, polarization, use_random, random_seed, jitter_amount)
    if ~isempty(random_seed), rng(random_seed); end
    
    if polarization < 0.05
        if use_random, theta_i = sort(rand(1, N) * 2 * pi);
        else, theta_i = linspace(0, 2*pi*(N-1)/N, N); end
    else
        power = 1 + polarization * 5;
        theta_fine = linspace(0, 2*pi, 5000);
        density = (sin(theta_fine).^2).^power + 1e-10;
        density = density / sum(density);
        cdf = cumsum(density); cdf = cdf / cdf(end);
        [cdf_u, idx] = unique(cdf);
        
        if use_random, u = sort(rand(1, N)); else, u = linspace(0, 1, N); end
        u = max(min(u, cdf_u(end)), cdf_u(1));
        
        theta_i = interp1(cdf_u, theta_fine(idx), u, 'linear', 'extrap');
        theta_i = mod(theta_i + jitter_amount * randn(size(theta_i)), 2*pi);
    end
end

function [delta, converged, n_iter, final_res] = solve_equilibrium(delta_init, r_theta, r0, h, k_i, L0, F0, bonded, F_inertial, max_iter, tol, params)
    delta = delta_init;
    damping = params.init_damp;
    converged = false;
    
    % Filter bonded fibers immediately for speed
    r_b = r_theta(bonded);
    k_b = k_i(bonded);
    L0_b = L0(bonded);
    
    if isempty(r_b)
        converged = true; n_iter = 0; final_res = 0; return;
    end
    
    for iter = 1:max_iter
        phi = atan((h + delta) ./ (r_b - r0));
        L_f = sqrt((r_b - r0).^2 + (h + delta).^2);
        F_i = F0 + k_b .* h .* (L_f ./ L0_b - 1);
        
        residual = sum(F_i .* sin(phi)) - F_inertial;
        
        if abs(residual) < tol
            converged = true; n_iter = iter; final_res = abs(residual); return;
        end
        
        % Analytical Derivatives
        dF = k_b .* h ./ L0_b .* ((h + delta) ./ L_f);
        dPhi = (r_b - r0) ./ ((r_b - r0).^2 + (h + delta).^2);
        dR = sum(dF .* sin(phi) + F_i .* cos(phi) .* dPhi);
        
        if abs(dR) < 1e-12, dR = sign(dR)*1e-12; if dR==0, dR=1e-12; end; end
        
        step = -residual / dR;
        
        % Adaptive Damping
        if params.use_adaptive
            delta_test = max(0, delta + damping * step);
            phi_t = atan((h + delta_test) ./ (r_b - r0));
            L_t = sqrt((r_b - r0).^2 + (h + delta_test).^2);
            F_t = F0 + k_b .* h .* (L_t ./ L0_b - 1);
            res_t = abs(sum(F_t .* sin(phi_t)) - F_inertial);
            
            if res_t < abs(residual)
                delta = delta_test;
                damping = min(damping * params.damp_inc, params.max_damp);
            else
                damping = max(damping * params.damp_dec, params.min_damp);
                % Conservative update option: delta = max(0, delta + damping * step);
            end
        else
            delta = max(0, delta + damping * step);
        end
    end
    n_iter = max_iter; final_res = abs(residual);
end

function [bonded_hist, F_hist, N_hist, S_hist, theta_i, r_theta, delta_hist, prob_debond_hist, conv_stats, k0_final] = ...
    simulate_full_history(a, b, r0, h, k0, F0, c, F_inertial, N, dt, n_steps, max_iter, tol, k_off0, N0, K_s, K_0, n_exp, m_exp, ...
    polarization, use_random, seed, jitter, check_conv, params)

    theta_i = generate_fiber_angles(N, polarization, use_random, seed, jitter);
    r_theta = a ./ sqrt(cos(theta_i).^2 + (a/b)^2 * sin(theta_i).^2);
    
    bonded_hist = false(n_steps, N);
    F_hist = zeros(n_steps, N); N_hist = zeros(n_steps, N); S_hist = zeros(n_steps, N);
    prob_debond_hist = zeros(n_steps, N); delta_hist = zeros(n_steps, 1);
    
    iter_hist = nan(n_steps, 1); res_hist = nan(n_steps, 1); conv_flags = false(n_steps, 1);
    
    % Auto-adjust loop
    k0_curr = k0;
    success = false;
    
    for attempt = 1:3
        L0 = sqrt((r_theta - r0).^2 + h^2);
        k_i = k0_curr * h ./ L0;
        bonded = true(1, N);
        delta = 0;
        
        fail = false;
        for t = 1:n_steps
            if sum(bonded) < 5, fail = false; break; end % Detached, but simulation valid
            
            [delta, conv, iter, res] = solve_equilibrium(delta, r_theta, r0, h, k_i, L0, F0, bonded, F_inertial, max_iter, tol, params);
            
            iter_hist(t) = iter; res_hist(t) = res; conv_flags(t) = conv;
            
            if ~conv && params.auto_adjust_k0 && attempt < 3
                fail = true; break;
            end
            
            % Update Forces
            phi = atan((h + delta) ./ (r_theta - r0));
            L_f = sqrt((r_theta - r0).^2 + (h + delta).^2);
            F_i = F0 + k_i .* h .* (L_f ./ L0 - 1);
            N_i = F_i .* sin(phi);
            S_i = F_i .* cos(phi);
            
            % Kinetics
            mask = bonded;
            k_off = zeros(1, N);
            if any(mask)
                term1 = (K_s^n_exp) ./ (K_s^n_exp + S_i(mask).^n_exp);
                term2 = 1 + c * (S_i(mask).^m_exp ./ (K_0^m_exp + S_i(mask).^m_exp));
                k_off(mask) = k_off0 * exp(N_i(mask)/N0 .* term1 .* term2);
            end
            prob = 1 - exp(-k_off * dt);
            
            bonded_hist(t, :) = bonded;
            F_hist(t, :) = F_i; N_hist(t, :) = N_i; S_hist(t, :) = S_i;
            prob_debond_hist(t, :) = prob;
            delta_hist(t) = delta;
            
            bonded(bonded & (rand(1, N) <= prob)) = false;
        end
        
        if fail
            fprintf('  Convergence fail at attempt %d. Reducing k0...\n', attempt);
            k0_curr = k0_curr * params.k0_reduction;
        else
            success = true; k0_final = k0_curr; break;
        end
    end
    
    conv_stats.iter_history = iter_hist;
    conv_stats.residual_history = res_hist;
    conv_stats.fail_rate = mean(~conv_flags);
    conv_stats.avg_iter = mean(iter_hist(~isnan(iter_hist)));
    conv_stats.avg_residual = mean(res_hist(~isnan(res_hist)));
end

function plot_timestep(ax, a, b, r0, theta_i, r_theta, bonded, prob, AR, t, dt)
    axes(ax); cla; hold on; axis equal;
    
    % Ellipse
    th = linspace(0, 2*pi, 200);
    plot(a*cos(th), b*sin(th), 'k-', 'LineWidth', 1.5);
    % Cap
    plot(r0*cos(th), r0*sin(th), 'k--', 'LineWidth', 1);
    
    % Fibers
    cmap = jet(256);
    if sum(bonded) > 0
        p_min = min(prob(bonded)); p_max = max(prob(bonded)) + 1e-10;
    else
        p_min = 0; p_max = 1;
    end
    caxis([p_min, p_max]); colormap(ax, cmap);
    
    for i = 1:length(theta_i)
        x0 = r0*cos(theta_i(i)); y0 = r0*sin(theta_i(i));
        x1 = r_theta(i)*cos(theta_i(i)); y1 = r_theta(i)*sin(theta_i(i));
        
        if bonded(i)
            val = (prob(i) - p_min) / (p_max - p_min);
            c_idx = round(val*255) + 1;
            plot([x0 x1], [y0 y1], 'Color', cmap(c_idx, :), 'LineWidth', 1.5);
        else
            plot([x0 x1], [y0 y1], 'Color', [0.8 0.8 0.8], 'LineStyle', ':');
        end
    end
    
    title(sprintf('AR=%.1f | t=%.2f | %d Bonded', AR, t*dt, sum(bonded)));
    xlim([-b*1.1 b*1.1]); ylim([-b*1.1 b*1.1]);
end