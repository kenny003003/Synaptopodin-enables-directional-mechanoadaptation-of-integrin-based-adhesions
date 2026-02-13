clear; close all; clc;

%% === PARAMETERS ===
% Conditions
C_VALUES = [1.3, 2.45];      % Control vs KO
LABELS = {'Control (C=1.3)', 'Synpo^{-/-} (C=2.45)'}; % Updated for LaTeX formatting

% Publication-Level Colors (RGB Triplets)
% Deep Blue for Control, Brick Red for KO
COLOR_MAP = {[179,229,252]/255,[251,154,154]/255}; 
FILL_MAP  = {[179,229,252]/255,[251,154,154]/255}; 

% Simulation Settings
N_TRIALS = 1000;             
N_FIBERS = 100;
MAX_TIMESTEPS = 250;
DT = 0.01;

% Physical Constants
F_BASE = 6;
F_RATIO = 3.5;               
F_INERTIAL = F_RATIO * F_BASE; 
AR_VAL = 1.4;                
r0 = 0.3; h = 0.1; k0 = 10000000; F0 = 0.0;
k_off0 = 0.8; N0 = 1; K_s = 1; K_0 = 1; n_exp = 1.02; m_exp = 2;

%% === SIMULATION LOOP ===
ar_results = cell(length(C_VALUES), 1);

if isempty(gcp('nocreate'))
    parpool; 
end

fprintf('Starting batch simulation (%d trials)...\n', N_TRIALS);
t_total = tic;

for c_idx = 1:length(C_VALUES)
    C_val = C_VALUES(c_idx);
    fprintf('Processing Condition %d...\n', c_idx);
    
    radius_base = 1.0;
    a = radius_base / sqrt(AR_VAL);
    b = radius_base * sqrt(AR_VAL); 
    
    ar_history_matrix = NaN(N_TRIALS, MAX_TIMESTEPS);
    
    parfor trial = 1:N_TRIALS
        trial_ar_history = NaN(1, MAX_TIMESTEPS);
        theta_i = generate_fiber_angles_uniform(N_FIBERS);
        
        r_theta = a ./ sqrt(cos(theta_i).^2 + (a/b)^2 * sin(theta_i).^2);
        L0 = sqrt((r_theta - r0).^2 + h^2);
        k_i = k0 * h ./ L0;
        
        bonded = true(1, N_FIBERS);
        delta = 0;
        
        for t = 1:MAX_TIMESTEPS
            if sum(bonded) < (N_FIBERS * 0.05)
                break; 
            end
            
            [delta, bonded] = solve_physics_step(delta, bonded, ...
                r_theta, r0, h, k_i, L0, F0, F_INERTIAL, theta_i, ...
                k_off0, N0, C_val, DT, K_s, K_0, n_exp, m_exp);
            
            if sum(bonded) >= 5
                 x_b = r_theta(bonded) .* cos(theta_i(bonded));
                 y_b = r_theta(bonded) .* sin(theta_i(bonded));
                 
                 pts = [x_b', y_b'];
                 G_tensor = cov(pts); 
                 eigenvalues = eig(G_tensor);
                 
                 if min(eigenvalues) > 1e-9
                     current_ar = sqrt(max(eigenvalues)) / sqrt(min(eigenvalues));
                 else
                     current_ar = NaN;
                 end
                 trial_ar_history(t) = current_ar;
            else
                trial_ar_history(t) = NaN;
            end
        end
        ar_history_matrix(trial, :) = trial_ar_history;
    end
    ar_results{c_idx} = ar_history_matrix;
end
toc(t_total);

%% === PUBLICATION PLOTTING ===
% Create figure with white background and specific aspect ratio
figure('Color', 'w', 'Position', [100 100 700 550]); 
hold on;

% Define Axes Properties upfront for clean look
ax = gca;
ax.FontSize = 16;
ax.LineWidth = 1.5;
ax.Box = 'on';
ax.TickDir = 'out';
ax.FontName = 'Arial'; % or 'Helvetica'
ax.XColor = [0.1 0.1 0.1];
ax.YColor = [0.1 0.1 0.1];

time_axis = (1:MAX_TIMESTEPS) * DT;

for c_idx = 1:length(C_VALUES)
    data = ar_results{c_idx}; 
    
    % Statistics
    mu = mean(data, 1, 'omitnan');
    sigma = std(data, 0, 1, 'omitnan');
    n_valid = sum(~isnan(data), 1);
    
    % 95% Confidence Interval
    sem = sigma ./ sqrt(n_valid);
    ci_upper = mu + 1.96 * sem;
    ci_lower = mu - 1.96 * sem;
    
    valid_mask = n_valid >= 10;
    
    t_plot = time_axis(valid_mask);
    upper_plot = ci_upper(valid_mask);
    lower_plot = ci_lower(valid_mask);
    mu_plot = mu(valid_mask);
    
    % 1. Draw Shaded Area (Confidence Interval)
    fill([t_plot, fliplr(t_plot)], [upper_plot, fliplr(lower_plot)], ...
         FILL_MAP{c_idx}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off'); % Turn off legend handle for fill
     
    % 2. Draw Mean Line
    plot(t_plot, mu_plot, 'Color', COLOR_MAP{c_idx}, 'LineWidth', 3, ...
         'DisplayName', LABELS{c_idx});

end

% Labels
xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Aspect Ratio', 'FontSize', 18, 'FontWeight', 'bold');

% Legend (Publication Style: No Box, specific location)
lgd = legend('Location', 'best');
lgd.EdgeColor = 'none'; % Remove box outline
lgd.FontSize = 14;

% Limits and Grid
ylim([1, 2.2]); 
xlim([0, MAX_TIMESTEPS*DT]);
% grid on; % Optional: Many high-impact journals prefer no grid, uncomment if needed

hold off;

%% === HELPER FUNCTIONS ===
function theta_i = generate_fiber_angles_uniform(N)
    theta_i = linspace(0, 2*pi, N+1);
    theta_i(end) = []; 
    theta_i = mod(theta_i + 0.001 * randn(1, N), 2*pi);
end

function [delta_new, bonded_new] = solve_physics_step(delta_start, bonded_mask, ...
    r_theta, r0, h, k_i, L0, F0, F_load, theta_i, ...
    k_off0, N0, C_val, dt, K_s, K_0, n_exp, m_exp)
    
    bonded_new = bonded_mask;
    delta_new = delta_start;
    damping = 0.3;
    
    r_b = r_theta(bonded_new);
    k_b = k_i(bonded_new);
    L0_b = L0(bonded_new);
    
    if isempty(r_b), return; end
    
    % Equilibrium Solver
    for i = 1:100
        phi = atan((h + delta_new) ./ (r_b - r0));
        L_f = sqrt((r_b - r0).^2 + (h + delta_new)^2);
        F_i = F0 + k_b .* h .* (L_f ./ L0_b - 1);
        
        residual = sum(F_i .* sin(phi)) - F_load;
        if abs(residual) < 1e-6, break; end
        
        dF = k_b .* h ./ L0_b .* ((h + delta_new) ./ L_f);
        dPhi = (r_b - r0) ./ ((r_b - r0).^2 + (h + delta_new).^2);
        dR = sum(dF .* sin(phi) + F_i .* cos(phi) .* dPhi);
        
        if abs(dR) < 1e-12, dR = sign(dR)*1e-12; end
        step = -residual / dR;
        delta_new = max(0, delta_new + damping * step);
    end
    
    % Kinetics
    phi_full = atan((h + delta_new) ./ (r_theta - r0));
    L_full = sqrt((r_theta - r0).^2 + (h + delta_new).^2);
    F_full = F0 + k_i .* h .* (L_full ./ L0 - 1);
    
    if any(bonded_new)
        N_force = F_full(bonded_new) .* sin(phi_full(bonded_new));
        S_force = F_full(bonded_new) .* cos(phi_full(bonded_new));
        
        term1 = (K_s^n_exp) ./ (K_s^n_exp + S_force.^n_exp);
        term2 = 1 + C_val .* (S_force.^m_exp ./ (K_0^m_exp + S_force.^m_exp));
        k_off = k_off0 .* exp(N_force./N0 .* term1 .* term2);
        
        prob = 1 - exp(-k_off * dt);
        broken_indices = find(bonded_new);
        is_broken = rand(size(broken_indices)) < prob;
        bonded_new(broken_indices(is_broken)) = false;
    end
end

%% === DATA EXPORT FOR PRISM ===
fprintf('Exporting data for GraphPad Prism...\n');

% Initialize variables for table storage
% Prism Format: Time, Mean(A), SEM(A), N(A), Mean(B), SEM(B), N(B)...
export_table = table();
export_table.Time_s = time_axis';

for c_idx = 1:length(C_VALUES)
    data = ar_results{c_idx}; 
    
    % Recalculate Statistics (Same as plotting section)
    mu = mean(data, 1, 'omitnan')';
    sigma = std(data, 0, 1, 'omitnan')';
    n_valid = sum(~isnan(data), 1)';
    sem = sigma ./ sqrt(n_valid);
    
    % Create safe column names
    col_base = regexprep(LABELS{c_idx}, '[^\w]', '_'); % Remove LaTeX/Special chars
    col_base = erase(col_base, {'__', '{', '}'});      % Cleanup
    
    % Add to table
    export_table.(['Mean_' col_base]) = mu;
    export_table.(['SEM_' col_base]) = sem;
    export_table.(['N_' col_base]) = n_valid;
end

% Write to Excel
filename = 'Prism_Export_AR_Data.xlsx';
writetable(export_table, filename);
fprintf('Success! Data saved to: %s\n', fullfile(pwd, filename));