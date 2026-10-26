clear; close all; clc;

%% === PARAMETERS ===
% Conditions
C_VALUES = [1.5, 2.7];      % Control vs KO
LABELS = {'Control (C=1.5)', 'Synpo^{-/-} (C=2.7)'}; % Updated for LaTeX formatting

% Publication-Level Colors (RGB Triplets)
% Deep Blue for Control, Brick Red for KO
COLOR_MAP = {[179,229,252]/255,[251,154,154]/255}; 
FILL_MAP  = {[179,229,252]/255,[251,154,154]/255}; 

% Simulation Settings
N_TRIALS = 400;             
N_FIBERS = 100;
MAX_TIMESTEPS = 250;
DT = 0.01;

% Physical Constants
F_BASE = 6;
F_RATIO = 3.5;               
F_INERTIAL = F_RATIO * F_BASE; 
AR_VAL = 1.4;                
r0 = 0.3; h = 0.1; k0 = 1000000; F0 = 2.5;
k_off0 = 0.8; N0 = 2.5; K_s = 2.5; K_0 = 2.5; n_exp = 1.1; m_exp = 4;

%% === SIMULATION LOOP ===
ar_results = cell(length(C_VALUES), 1);

fprintf('Starting batch simulation (%d trials)...\n', N_TRIALS);
t_total = tic;

for c_idx = 1:length(C_VALUES)
    C_val = C_VALUES(c_idx);
    fprintf('Processing Condition %d...\n', c_idx);
    
    radius_base = 1.0;
    a = radius_base / sqrt(AR_VAL);
    b = radius_base * sqrt(AR_VAL); 
    
    ar_history_matrix = NaN(N_TRIALS, MAX_TIMESTEPS);
    
    for trial = 1:N_TRIALS
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
figure('Color','w','Position',[100 100 780 460],'Visible','off'); 
hold on;

% Define Axes Properties upfront for clean look
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.0;
ax.Box = 'on';
ax.TickDir = 'in';
ax.FontName = 'Arial'; % or 'Helvetica'
ax.XColor = [0.1 0.1 0.1];
ax.YColor = [0.1 0.1 0.1];

time_axis = (1:MAX_TIMESTEPS) / MAX_TIMESTEPS;

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
    plot(t_plot, mu_plot, 'Color', COLOR_MAP{c_idx}, 'LineWidth', 2.2, ...
         'DisplayName', LABELS{c_idx});

end

% Labels
%% xlabel removed
%% ylabel removed

% Legend (Publication Style: No Box, specific location)
lgd = legend('Location', 'northwest');
lgd.EdgeColor = 'none'; % Remove box outline
lgd.FontSize = 12;

% Limits and Grid
ylim([1.2 2]); yticks([1.2 1.4 1.6 1.8 2]); 
xlim([0 1]); set(ax,'XTickLabel',[],'YTickLabel',[]);
% grid on; % Optional: Many high-impact journals prefer no grid, uncomment if needed

hold off;

exportgraphics(gcf, "aspect_ratio_sweep.png", "Resolution", 170);
fprintf("Saved aspect_ratio_sweep.png\n");


%% Export two-curve data to Excel
export_T = table();
export_T.Normalized_time = time_axis(:);
condtags = {'Control_C1p5','Synpo_C2p7'};
for ci = 1:length(C_VALUES)
    dd = ar_results{ci};
    mu = mean(dd,1,"omitnan").';
    sg = std(dd,0,1,"omitnan").';
    nv = sum(~isnan(dd),1).';
    sem = sg./sqrt(nv);
    tg = condtags{ci};
    export_T.([tg '_mean']) = mu;
    export_T.([tg '_SEM']) = sem;
    export_T.([tg '_CI95low']) = mu - 1.96*sem;
    export_T.([tg '_CI95high']) = mu + 1.96*sem;
end
writetable(export_T, "panelB_curve_data.xlsx");
fprintf("Saved panelB_curve_data.xlsx (%d rows)\n", height(export_T));


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
