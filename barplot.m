clear; close all; clc;

%% === PARAMETERS (Matched to parameter_sweep_v3) ===
% 1. Force Groups
F_BASE = 5;
FORCE_RATIOS = [1, 3.5, 15]; 
FORCE_LABELS = {'Low Force', 'Med Force', 'High Force'};

% 2. Biological Conditions 
% Values taken from parameter_sweep_v3: C_VALUES=[1.3, 2.45], AR_VALUES=[1.4]
% Structure: [C_value, Aspect_Ratio]
CONDITIONS = [
    1.3,  1.4;   % Condition 1 (Control)
    2.45, 1.4    % Condition 2 (Synops KO)
];
COND_NAMES = {'Control (C=1.3)', 'Synops KO (C=2.45)'};
COND_COLORS = [0.4 0.7 1.0;  % Light Blue 
               1.0 0.4 0.4]; % Light Red 

% 3. Simulation Control
N_SETS = 20;                % Sets of data points (for error bars)
N_CELLS_PER_SET = 50;       % Cells per set (Total = 20*50 = 1000 cells)
N_FIBERS = 50;              % *** FIXED: Changed from 100 to 50 to match v3 ***
DETACH_THRESHOLD = 0.05;
MAX_TIMESTEPS = 250;
DT = 0.01;

% 4. Physical Constants (Matched to v3)
r0 = 0.3; h = 0.1; k0 = 10000000; F0 = 0.0;
k_off0 = 0.8; N0 = 1; K_s = 1; K_0 = 1; n_exp = 1.05; m_exp = 2;

%% === SIMULATION LOOP ===
% Storage: attachment_data(Force_Idx, Condition_Idx, Set_Idx)
attachment_data = zeros(length(FORCE_RATIOS), size(CONDITIONS, 1), N_SETS);

% Check for parallel pool
if isempty(gcp('nocreate'))
    parpool; 
end

fprintf('=== Starting Simulation (Attachment Rate) ===\n');
fprintf('Running %d sets of %d cells for each condition...\n', N_SETS, N_CELLS_PER_SET);
total_start = tic;

for f_idx = 1:length(FORCE_RATIOS)
    F_ratio = FORCE_RATIOS(f_idx);
    F_inertial = F_ratio * F_BASE;
    
    fprintf('\nProcessing Force Ratio: %.1f ...\n', F_ratio);
    
    for c_idx = 1:size(CONDITIONS, 1)
        % Extract Condition Parameters
        C_val = CONDITIONS(c_idx, 1);
        AR_val = CONDITIONS(c_idx, 2);
        
        % Geometry Setup (Constant Area)
        % Logic from v3: AR > 1.5 is polarized
        if AR_val > 1.5
            polarization = 0.7; 
        else
            polarization = 0.0; 
        end
        
        radius_base = 1.0;
        a_fixed = radius_base / sqrt(AR_val);
        b_axis = a_fixed * AR_val;
        
        % Temporary storage for parallel results
        results_temp = zeros(1, N_SETS);
        
        % Run sets in Parallel
        parfor set_k = 1:N_SETS
            % Run N cells for this set
            detached_count = 0;
            
            for cell_k = 1:N_CELLS_PER_SET
                % 1. Generate Fibers
                % (RNG is handled automatically by parfor)
                theta_i = generate_fiber_angles(N_FIBERS, polarization);
                
                % Geometry
                r_theta = a_fixed ./ sqrt(cos(theta_i).^2 + (a_fixed/b_axis)^2 * sin(theta_i).^2);
                L0 = sqrt((r_theta - r0).^2 + h^2);
                k_i = k0 * h ./ L0;
                
                % 2. Initialize
                bonded = true(1, N_FIBERS);
                delta = 0;
                is_detached = false;
                
                % 3. Time Steps
                for t = 1:MAX_TIMESTEPS
                    if sum(bonded) < (N_FIBERS * DETACH_THRESHOLD)
                        is_detached = true;
                        break;
                    end
                    
                    % Solver
                    [delta, ~] = solve_equilibrium(delta, r_theta, r0, h, k_i, L0, F0, bonded, F_inertial);
                    
                    % Forces & Kinetics
                    phi_i = atan((h + delta) ./ (r_theta - r0));
                    L_f = sqrt((r_theta - r0).^2 + (h + delta).^2);
                    F_i = F0 + k_i .* h .* (L_f ./ L0 - 1);
                    N_i = F_i .* sin(phi_i);
                    S_i = F_i .* cos(phi_i);
                    
                    k_off = zeros(1, N_FIBERS);
                    mask = bonded;
                    if any(mask)
                        term1 = (K_s^n_exp) ./ (K_s^n_exp + S_i(mask).^n_exp);
                        term2 = 1 + C_val .* (S_i(mask).^m_exp ./ (K_0^m_exp + S_i(mask).^m_exp));
                        k_off(mask) = k_off0 .* exp(N_i(mask)./N0 .* term1 .* term2);
                    end
                    
                    prob = 1 - exp(-k_off * DT);
                    bonded(bonded & (rand(1, N_FIBERS) <= prob)) = false;
                end
                
                if is_detached
                    detached_count = detached_count + 1;
                end
            end
            
            % Store ATTACHMENT rate (%) for this set
            results_temp(set_k) = ((N_CELLS_PER_SET - detached_count) / N_CELLS_PER_SET) * 100; 
        end
        
        % Transfer temp results to main storage
        attachment_data(f_idx, c_idx, :) = results_temp;
        fprintf('   %s: Mean Attachment = %.1f%%\n', COND_NAMES{c_idx}, mean(results_temp));
    end
end
toc(total_start);

%% === VISUALIZATION ===
figure('Position', [100 100 1000 600]);

% Calculate Statistics
means = mean(attachment_data, 3); % Mean across sets
stds = std(attachment_data, 0, 3); % STD across sets

% 1. Draw Bar Plot
b = bar(means, 'grouped');
b(1).FaceColor = COND_COLORS(1,:);
b(2).FaceColor = COND_COLORS(2,:);
b(1).DisplayName = COND_NAMES{1};
b(2).DisplayName = COND_NAMES{2};
hold on;

% 2. Add Error Bars and Scatter Points
[num_groups, num_bars] = size(means);
x_centers = zeros(num_groups, num_bars);

for i = 1:num_bars
    % Get X coordinates of the bars
    x_centers(:,i) = b(i).XEndPoints;
    
    % Plot Error Bars
    errorbar(x_centers(:,i), means(:,i), stds(:,i), 'k.', 'LineWidth', 1.5, 'CapSize', 10);
    
    % Plot Scatter Points (Jittered)
    for j = 1:num_groups
        % Extract data for this specific bar
        data_points = squeeze(attachment_data(j, i, :));
        
        % Jitter X coordinates
        x_jitter = x_centers(j,i) + (rand(size(data_points))-0.5) * 0.15;
        
        % Scatter
        scatter(x_jitter, data_points, 20, 'k', 'o', 'filled', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

% 3. Styling
ylabel('Attachment Rate (%)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Force Conditions', 'FontSize', 14, 'FontWeight', 'bold');
title('Cell Attachment: Control vs Synops KO', 'FontSize', 16);
xticks(1:length(FORCE_RATIOS));
xticklabels(FORCE_LABELS);
ylim([0 105]); 
legend('Location', 'northeast', 'FontSize', 12);
grid on;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);
hold off;


%% === HELPER FUNCTIONS ===

function theta_i = generate_fiber_angles(N, polarization)
    if polarization < 0.05
        theta_i = sort(rand(1, N) * 2 * pi);
    else
        power = 1 + polarization * 5;
        theta_fine = linspace(0, 2*pi, 1000);
        density = (sin(theta_fine).^2).^power + 1e-10;
        density = density / sum(density);
        cdf = cumsum(density); cdf = cdf / cdf(end);
        [cdf_u, idx] = unique(cdf);
        u = sort(rand(1, N));
        u = max(min(u, cdf_u(end)), cdf_u(1));
        theta_i = interp1(cdf_u, theta_fine(idx), u, 'linear');
        theta_i = mod(theta_i + 0.005 * randn(1, N), 2*pi);
    end
end

function [delta, converged] = solve_equilibrium(delta_init, r_theta, r0, h, k_i, L0, F0, bonded, F_inertial)
    delta = delta_init; damping = 0.3; max_iter = 100; tol = 1e-6; converged = false;
    r_b = r_theta(bonded); k_b = k_i(bonded); L0_b = L0(bonded);
    if isempty(r_b), converged=true; return; end
    for i = 1:max_iter
        phi = atan((h+delta)./(r_b-r0)); Lf = sqrt((r_b-r0).^2+(h+delta).^2);
        Fi = F0 + k_b.*h.*(Lf./L0_b-1);
        residual = sum(Fi.*sin(phi)) - F_inertial;
        if abs(residual)<tol, converged=true; return; end
        dF = k_b.*h./L0_b.*((h+delta)./Lf); dPhi = (r_b-r0)./((r_b-r0).^2+(h+delta).^2);
        dR = sum(dF.*sin(phi) + Fi.*cos(phi).*dPhi);
        if abs(dR)<1e-12, dR=sign(dR)*1e-12; if dR==0, dR=1e-12; end; end
        delta = max(0, delta - damping*(residual/dR));
    end
end