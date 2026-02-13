%% === COMPLETE SIMULATION & BUTTERFLY VISUALIZATION (With SEM) ===
clear; close all; clc;

% --- 1. SETUP & CONSTANTS ---
% Conditions
C_VALUES = [1.3, 2.45];
% Colors for Control (Left)
COLOR_CTL_FILL = [0.6 0.6 0.6]; % Grey
COLOR_CTL_LINE = [0 0 0];       % Black
% Colors for Synpo-/- (Right)
COLOR_KO_FILL  = [0.6 0.8 0.6]; % Light Green
COLOR_KO_LINE  = [0 0.5 0];     % Dark Green

% Simulation Parameters
N_TRIALS = 100;               
N_FIBERS = 100;
MAX_TIMESTEPS = 250;
DT = 0.01;
DETACH_THRESHOLD = 0.05; 

% Physical Constants
F_BASE = 6;
F_RATIO = 3.5;
F_INERTIAL = F_RATIO * F_BASE; 
AR_VAL = 1.4;
r0 = 0.3; h = 0.1; k0 = 10000000; F0 = 0.0;
k_off0 = 0.8; N0 = 1; K_s = 1; K_0 = 1; n_exp = 1.05; m_exp = 2;

% Rose Plot Targets (0 to 180)
target_angles_deg = [0, 30, 60, 90, 120, 150, 180];
target_angles_rad = deg2rad(target_angles_deg);
bin_width_rad = deg2rad(15);

% Geometry for Ellipse
radius_base = 1.0;
a = radius_base * sqrt(AR_VAL);  
b = radius_base / sqrt(AR_VAL);  

% Storage
all_rose_data = cell(1, 2);

%% --- 2. SIMULATION LOOP ---
if isempty(gcp('nocreate')), parpool; end
for c_idx = 1:length(C_VALUES)
    C_CURRENT = C_VALUES(c_idx);
    fprintf('Running %d trials for C=%.2f...\n', N_TRIALS, C_CURRENT);
    
    current_data = NaN(length(target_angles_deg), N_TRIALS);
    
    parfor trial = 1:N_TRIALS
        % Initialize
        theta_i = linspace(0, 2*pi, N_FIBERS+1);
        theta_i(end) = []; 
        theta_i = mod(theta_i + 0.001 * randn(1, N_FIBERS), 2*pi);
        
        r_theta = a ./ sqrt(cos(theta_i).^2 + (a/b)^2 * sin(theta_i).^2);
        L0 = sqrt((r_theta - r0).^2 + h^2);
        k_i = k0 * h ./ L0;
        
        bonded = true(1, N_FIBERS);
        delta = 0;
        is_detached = false;
        
        % Physics
        for t = 1:MAX_TIMESTEPS
            if sum(bonded) < (N_FIBERS * DETACH_THRESHOLD)
                is_detached = true; break; 
            end
            
            r_b = r_theta(bonded); k_b = k_i(bonded); L0_b = L0(bonded);
            if ~isempty(r_b)
                for iter = 1:100
                    phi = atan((h + delta) ./ (r_b - r0));
                    L_f = sqrt((r_b - r0).^2 + (h + delta).^2);
                    F_i = F0 + k_b .* h .* (L_f ./ L0_b - 1);
                    residual = sum(F_i .* sin(phi)) - F_INERTIAL;
                    if abs(residual) < 1e-6, break; end
                    dF = k_b .* h ./ L0_b .* ((h + delta) ./ L_f);
                    dPhi = (r_b - r0) ./ ((r_b - r0).^2 + (h + delta).^2);
                    dR = sum(dF .* sin(phi) + F_i .* cos(phi) .* dPhi);
                    if abs(dR) < 1e-12, dR = sign(dR)*1e-12; end
                    delta = max(0, delta - 0.3 * (residual / dR));
                end
            end
            
            % Kinetics
            phi_full = atan((h + delta) ./ (r_theta - r0));
            L_full = sqrt((r_theta - r0).^2 + (h + delta).^2);
            F_full = F0 + k_i .* h .* (L_full ./ L0 - 1);
            
            if any(bonded)
                N_f = F_full(bonded) .* sin(phi_full(bonded));
                S_f = F_full(bonded) .* cos(phi_full(bonded));
                term1 = (K_s^n_exp) ./ (K_s^n_exp + S_f.^n_exp);
                term2 = 1 + C_CURRENT .* (S_f.^m_exp ./ (K_0^m_exp + S_f.^m_exp));
                k_off = k_off0 .* exp(N_f./N0 .* term1 .* term2);
                prob = 1 - exp(-k_off * DT);
                
                idx_b = find(bonded);
                broken = rand(size(idx_b)) < prob;
                bonded(idx_b(broken)) = false;
            end
        end
        
        % Data Collection
        if ~is_detached
            all_thetas_folded = mod(theta_i, pi);
            binned_r = NaN(1, 7);
            for k = 1:7
                t_target = target_angles_rad(k);
                dist = abs(all_thetas_folded - t_target);
                bin_mask_all = dist <= bin_width_rad;
                if any(bin_mask_all)
                    total_count = sum(bin_mask_all);
                    surviving_mask = bin_mask_all & bonded;
                    sum_surviving_lengths = sum(r_theta(surviving_mask));
                    binned_r(k) = sum_surviving_lengths / total_count;
                else
                    binned_r(k) = NaN; 
                end
            end
            current_data(:, trial) = binned_r';
        end
    end
    all_rose_data{c_idx} = current_data;
end

%% === 3. VISUALIZATION (Clean Wedges + SEM Bars) ===
figure('Color', 'w', 'Position', [50 100 1200 600]);

% Error Bar Params
capWidth = deg2rad(2); % Width of the caps on error bars

% Helper to get theoretical radius
get_r_ellipse = @(ang) a ./ sqrt(cos(ang).^2 + (a/b).^2 .* sin(ang).^2);

for c_idx = 1:2
    % 1. Setup Subplot
    subplot(1, 2, c_idx);
    hold on; axis equal; axis off;
    
    % --- LAYER 1: EQUATOR LINE (Background) ---
    plot([-1.6 1.6], [0 0], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    
    % 2. Determine Colors
    if c_idx == 1
        % === CONTROL CASE ===
        titleStr = 'Control';
        topFill  = COLOR_CTL_FILL; topLine  = COLOR_CTL_LINE; 
        botFill  = COLOR_KO_FILL;  botLine  = COLOR_KO_LINE;   
    else
        % === SYNPO-/- CASE ===
        titleStr = 'Synpo^{-/-}';
        topFill  = COLOR_CTL_FILL; topLine  = COLOR_CTL_LINE;   
        botFill  = COLOR_KO_FILL;  botLine  = COLOR_KO_LINE;   
    end
    
    rose_data = all_rose_data{c_idx};
    
    % 3. Loop Through Angles (0 to 180)
    for i = 1:7
        theta_target = target_angles_rad(i);
        
        %% --- PART 1: TOP (INITIAL STATE) ---
        if i == 1 % 0 degrees
            t_start_top = 0; t_end_top = bin_width_rad;
        elseif i == 7 % 180 degrees
            t_start_top = pi - bin_width_rad; t_end_top = pi;
        else 
            t_start_top = theta_target - bin_width_rad; 
            t_end_top = theta_target + bin_width_rad;
        end
        sample_angles_top = linspace(t_start_top, t_end_top, 50);
        r_avg_top = mean(get_r_ellipse(sample_angles_top))/7.5; 
        
        [twx, twy] = pol2cart([t_start_top, t_end_top, 0], [r_avg_top, r_avg_top, 0]);
        patch(twx, twy, topFill, 'EdgeColor', topLine, ...
              'FaceAlpha', 0.8, 'LineWidth', 0.5);
        
        %% --- PART 2: BOTTOM (FORCED STATE - Wedges + SEM) ---
        theta_bottom = -theta_target; 
        rData = rose_data(i, :);
        rData = rData(~isnan(rData)); 
        
        if ~isempty(rData)
            % Stats
            meanVal = mean(rData);
            nVal = length(rData);
            semVal = std(rData) / sqrt(nVal); % Standard Error
            
            % A. Define Wedge Boundaries
            if i == 1 % Bottom 0 deg
                 t_start_bot = -bin_width_rad; t_end_bot = 0;
            elseif i == 7 % Bottom 180 deg
                 t_start_bot = -pi; t_end_bot = -pi + bin_width_rad;
            else
                 t_start_bot = theta_bottom - bin_width_rad; 
                 t_end_bot = theta_bottom + bin_width_rad;
            end
            
            % Draw Filled Wedge
            [bwx, bwy] = pol2cart([t_start_bot, t_end_bot, 0], [meanVal, meanVal, 0]);
            patch(bwx, bwy, botFill, 'EdgeColor', botLine, ...
                  'FaceAlpha', 0.8, 'LineWidth', 0.5);
            
            % B. Draw SEM Bars (Radial)
            lowerSE = meanVal - semVal;
            upperSE = meanVal + semVal;
            
            % Angle for the bar (center of the standard wedge)
            t_bar = theta_bottom; 
            
            % Main radial line
            [lx1, ly1] = pol2cart(t_bar, lowerSE);
            [lx2, ly2] = pol2cart(t_bar, upperSE);
            plot([lx1, lx2], [ly1, ly2], '-', 'Color', botLine, 'LineWidth', 1.5);
            
            % Caps (Tangential lines)
            % Cap at Lower SE
            [cx1a, cy1a] = pol2cart(t_bar - capWidth, lowerSE);
            [cx1b, cy1b] = pol2cart(t_bar + capWidth, lowerSE);
            plot([cx1a, cx1b], [cy1a, cy1b], '-', 'Color', botLine, 'LineWidth', 1.5);
            
            % Cap at Upper SE
            [cx2a, cy2a] = pol2cart(t_bar - capWidth, upperSE);
            [cx2b, cy2b] = pol2cart(t_bar + capWidth, upperSE);
            plot([cx2a, cx2b], [cy2a, cy2b], '-', 'Color', botLine, 'LineWidth', 1.5);
        end
    end
    
    % Add Center Hub
    plot(0,0, '.', 'Color', botLine, 'MarkerSize', 20);
    
    % Limits and Title
    title(titleStr, 'FontSize', 14, 'FontWeight', 'bold');
    xlim([-1.6, 1.6]); 
    ylim([-1.6, 1.6]); 
end