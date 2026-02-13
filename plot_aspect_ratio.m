%% === PUBLICATION PLOTTING ===
% Create figure with white background and specific aspect ratio
figure('Color', 'w', 'Position', [100 100 700 550]); 
hold on;

% --- AXIS CONFIGURATION ---
ax = gca;
ax.FontSize = 16;
ax.LineWidth = 1.5;
ax.Box = 'on';
ax.TickDir = 'out';
ax.FontName = 'Arial'; 
ax.XColor = [0.1 0.1 0.1];
ax.YColor = [0.1 0.1 0.1];
ylim=([1.1,2.1]);
set(gca, 'TickDir', 'in')
% --- CUSTOM TICK SETTINGS (EDIT HERE) ---
% Define the exact locations and text labels you want
x_tick_vals = [0, 0.5, 1.0, 1.5, 2.0, 2.5]; 
x_tick_labs = {'', '', '', '', '', ''};

y_tick_vals = [1.2, 1.4, 1.6, 1.8, 2.0];
y_tick_labs = {'1.2', '', '', '', '2'};

% Apply the ticks
xticks(x_tick_vals);
xticklabels(x_tick_labs);
yticks(y_tick_vals);
yticklabels(y_tick_labs);

% --- PLOTTING DATA ---
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
         'HandleVisibility', 'off'); 
     
    % 2. Draw Mean Line
    plot(t_plot, mu_plot, 'Color', COLOR_MAP{c_idx}, 'LineWidth', 3, ...
         'DisplayName', LABELS{c_idx});
end



% Legend
lgd = legend('Location', 'best');
lgd.EdgeColor = 'none'; 
lgd.FontSize = 14;

% Limits (Ensure limits cover your custom ticks)
xlim([min(x_tick_vals), max(x_tick_vals)]);

hold off;