close all, clear all, clc;

%% Two-Cable Triangle Configuration (Matched to Parameter Sweep)
% Represents 2 symmetric fibers from the full cell model
% Direct calculation - no iterative methods

%% Parameters (Matched to parameter_sweep_v2.m)
% Cell Geometry
r_cap = 0.3;             % Radius of actin cap (r0 in sweep)
r_cell = 1.0;            % Radius of cell boundary (r_theta for AR=1)
h = 0.05;                % Height of actin cap

% The "Triangle" base width is the horizontal span of the fiber
r_base = r_cell - r_cap; % Span from cap edge to cell edge (0.7)

% Fiber stiffness
k = 1000;            % Stiffness (k0 in sweep)

% Off-rate parameters (from sweep)
c = 1.0;                 % Off-rate parameter (C_VALUES)
K_0 = 1;                 % Catch-bond parameter (Km)
K_s = 1;                 % Shear force scale (Ks)
N0 = 1;                  % Characteristic normal force
k_off0 = 1;              % Base off-rate
n_exp = 1.5;               % Exponent for shear term
m_exp = 2;               % Exponent for catch-bond term
dt = 0.01;               % Time step

% Scaling for comparison
N_fibers_cell = 100;     % Number of fibers in full model
N_fibers_model = 2;      % Number of fibers in this model

% Initial cable length
L0 = sqrt(r_base^2 + h^2);

fprintf('=== Two-Cable Triangle Analysis (Matched to Sweep) ===\n');
fprintf('Geometry:\n');
fprintf('  Cap Radius r0 = %.2f\n', r_cap);
fprintf('  Cell Radius r_theta = %.2f\n', r_cell);
fprintf('  Fiber Horizontal Span = %.2f\n', r_base);
fprintf('  Height h = %.2f\n', h);
fprintf('  Stiffness k = %.1e\n', k);
fprintf('  Initial Angle = %.2f degrees\n\n', atan(h/r_base)*180/pi);

%% Direct calculation for range of displacements
% Sweep displacement delta
delta = linspace(0.022, 0.05, 1000);  

% Pre-allocate
L = zeros(size(delta));           
theta = zeros(size(delta));       
T = zeros(size(delta));           
F_applied = zeros(size(delta));   % Force on the 2-cable system
F_cell_equiv = zeros(size(delta));% Equivalent force on 100-fiber cell
N = zeros(size(delta));           
S = zeros(size(delta));           
k_off = zeros(size(delta));       
P_fail = zeros(size(delta));      

fprintf('Calculating force balance...\n');

for i = 1:length(delta)
    d = delta(i);
    
    % Current geometry (using r_base as the horizontal run)
    L(i) = sqrt(r_base^2 + (h + d)^2);       
    theta(i) = atan((h + d) / r_base);       
    
    % Tension T = k * (L - L0)
    T(i) = k * (L(i) - L0);                  
    
    % Vertical Equilibrium: F = 2 * T * sin(theta)
    F_applied(i) = 2 * T(i) * sin(theta(i));
    
    % Equivalent force if this were the full 100-fiber cell
    % Scale by ratio of fiber counts (x50)
    F_cell_equiv(i) = F_applied(i) * (N_fibers_cell / N_fibers_model);
    
    % Components per fiber
    N(i) = T(i) * sin(theta(i));             
    S(i) = T(i) * cos(theta(i));             
    
    % Off-rate calculation (per fiber)
    if N(i) > 0
        term1 = (K_s^n_exp) / (K_s^n_exp + S(i)^n_exp);
        term2 = 1 + c * (S(i)^m_exp / (K_0^m_exp + S(i)^m_exp));
        k_off(i) = k_off0 * exp((N(i)/N0) * term1 * term2);
    else
        k_off(i) = k_off0;
    end
    
    % Debonding probability
    P_fail(i) = 1 - exp(-k_off(i) * dt);
end

fprintf('Calculation complete!\n\n');

%% Comprehensive Plotting
figure('Position', [150 150 1400 900]);

% Plot 1: Displacement vs Force
subplot(2,3,1)
plot(F_cell_equiv, delta, 'k-', 'LineWidth', 2.5)
xlabel('Equivalent Cell Force F_{total}', 'FontSize', 12)
ylabel('Displacement \delta', 'FontSize', 12)
title('Displacement Response', 'FontWeight', 'bold', 'FontSize', 13)
grid on
set(gca, 'LineWidth', 1.5, 'FontSize', 11)

% Plot 2: Angle vs Force
subplot(2,3,2)
plot(F_cell_equiv, theta, 'k-', 'LineWidth', 2.5)
xlabel('Equivalent Cell Force F_{total}', 'FontSize', 12)
ylabel('Cable Angle \theta (degrees)', 'FontSize', 12)
title('Cable Angle Evolution', 'FontWeight', 'bold', 'FontSize', 13)
grid on
set(gca, 'LineWidth', 1.5, 'FontSize', 11)

% Plot 3: N/S Ratio vs Force
subplot(2,3,3)
plot(F_cell_equiv, N./S, 'k-', 'LineWidth', 2.5)
xlabel('Equivalent Cell Force F_{total}', 'FontSize', 12)
ylabel('N/S Ratio', 'FontSize', 12)
title('Normal/Shear Force Ratio', 'FontWeight', 'bold', 'FontSize', 13)
grid on
set(gca, 'LineWidth', 1.5, 'FontSize', 11)

% Plot 4: Force Components vs Applied Force
subplot(2,3,4)
plot(F_cell_equiv, N, 'r-', 'LineWidth', 2.5)
hold on
plot(F_cell_equiv, S, 'b-', 'LineWidth', 2.5)
plot(F_cell_equiv, T, 'k--', 'LineWidth', 2)
xlabel('Equivalent Cell Force F_{total}', 'FontSize', 12)
ylabel('Force (per cable)', 'FontSize', 12)
legend('Normal (N)', 'Shear (S)', 'Tension (T)', 'Location', 'northwest', 'FontSize', 10)
title('Force Components per Fiber', 'FontWeight', 'bold', 'FontSize', 13)
grid on
set(gca, 'LineWidth', 1.5, 'FontSize', 11)

% Plot 5: Off-rate vs Force
subplot(2,3,5)
semilogy(F_cell_equiv, k_off, 'k-', 'LineWidth', 2.5)
xlabel('Equivalent Cell Force F_{total}', 'FontSize', 12)
ylabel('k_{off}', 'FontSize', 12)
title('Off-Rate (log scale)', 'FontWeight', 'bold', 'FontSize', 13)
grid on
set(gca, 'LineWidth', 1.5, 'FontSize', 11)

% Plot 6: Failure Probability vs Force
subplot(2,3,6)
plot(F_cell_equiv, P_fail, 'k-', 'LineWidth', 2.5)
xlabel('Equivalent Cell Force F_{total}', 'FontSize', 12)
ylabel(sprintf('P(failure) in \\Deltat=%.2fs', dt), 'FontSize', 12)
title('Debonding Probability', 'FontWeight', 'bold', 'FontSize', 13)
grid on
set(gca, 'LineWidth', 1.5, 'FontSize', 11)

sgtitle(sprintf('Analytical Model (Matched to Sweep AR=1, C=%.1f)', c), ...
    'FontSize', 16, 'FontWeight', 'bold')