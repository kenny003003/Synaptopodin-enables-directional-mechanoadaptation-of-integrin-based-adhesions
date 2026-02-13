% Prompt user to select an Excel file
[file, path] = uigetfile('*.xlsx', 'Select the Excel File');
if isequal(file, 0)
    disp('User canceled the file selection');
    return;
end

% Read the first 7 rows of data
filename = fullfile(path, file);
data = readmatrix(filename, 'Range', 'A1:Z7');  % Adjust range as needed

% Prompt user for color choice
colorChoice = questdlg('Choose the color for the plot:', ...
    'Color Choice', ...
     'Blue', 'Red','Green', 'Black');  % Default: Black

% Set the chosen color
switch colorChoice
    
    case 'Red'
        plotColor = [1, 0, 0];
    case 'Blue'
        plotColor = [0, 0, 1];
    case 'Green'
        plotColor = [0, 0.5, 0];
    otherwise
        plotColor = [0, 0, 0];
end

% Ask whether to use median or mean for the central line
centerStatChoice = questdlg('Use median or mean for the central line?', ...
    'Central Tendency', ...
    'Median', 'Mean', 'Median');
useMedian = strcmp(centerStatChoice, 'Median');

% Ask whether to show a trend line
showTrendLine = questdlg('Do you want to show a trend line?', ...
    'Trend Line', ...
    'Yes', 'No', 'Yes');
useTrendLine = strcmp(showTrendLine, 'Yes');

% Define angles in degrees and convert to radians
angles_deg = [0, 30, 60, 90, 120, 150, 180];
angles_rad = deg2rad(angles_deg);

% Prepare center values for trend line
centerVals = zeros(1, 7);

% Create polar axes
figure;
pax = polaraxes;
hold on;

% Configure polar axes
pax.ThetaLim = [0 180];
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'bottom';
pax.RTick = [];
pax.ThetaColor = [0.6 0.6 0.6];
pax.RColor = 'none';
pax.ThetaTickLabel = {};

% Boxplot styling
boxWidth = deg2rad(4);
capWidth = boxWidth / 2;

for i = 1:7
    rData = data(i, :);
    theta = angles_rad(i);

    % Scatter raw data points with jitter in angle
    jitter = (rand(1, length(rData)) - 0.5) * boxWidth * 1.5;
    polarplot(theta + jitter, rData, 'o', ...
        'MarkerSize', 1, ...
        'MarkerEdgeColor', plotColor, ...
        'MarkerFaceColor', plotColor, ...
        'LineWidth', 0.5);

    % Boxplot statistics
    q1 = prctile(rData, 28);
    q3 = prctile(rData, 72);
    iqr = q3 - q1;
    lowerWhisker = max(min(rData), q1 - 1.1 * iqr);
    upperWhisker = min(max(rData), q3 + 1.1 * iqr);

    % Central line value
    if useMedian
        centerVal = prctile(rData, 50);
    else
        centerVal = mean(rData, 'omitnan');
    end
    centerVals(i) = centerVal;

    % Draw box
    thetaBox = [theta - boxWidth, theta + boxWidth, theta + boxWidth, theta - boxWidth, theta - boxWidth];
    rBox = [q1, q1, q3, q3, q1];
    polarplot(thetaBox, rBox, 'Color', plotColor, 'LineWidth', 1.5);

    % Central line (median or mean)
    polarplot([theta - boxWidth, theta + boxWidth], [centerVal, centerVal], ...
        'Color', plotColor, 'LineWidth', 3);

    % Whiskers
    polarplot([theta, theta], [lowerWhisker, q1], 'Color', plotColor, 'LineStyle', '--', 'LineWidth', 1);
    polarplot([theta, theta], [q3, upperWhisker], 'Color', plotColor, 'LineStyle', '--', 'LineWidth', 1);

    % Whisker caps
    polarplot([theta - capWidth, theta + capWidth], [lowerWhisker, lowerWhisker], 'Color', plotColor, 'LineWidth', 1);
    polarplot([theta - capWidth, theta + capWidth], [upperWhisker, upperWhisker], 'Color', plotColor, 'LineWidth', 1);
end

% ---- Curved trend line using interpolation ----
if useTrendLine
    thetaFine = linspace(min(angles_rad), max(angles_rad), 200);
    centerValsSmooth = interp1(angles_rad, centerVals, thetaFine, 'pchip');
    polarplot(thetaFine, centerValsSmooth, '-', 'Color', plotColor, 'LineWidth', 2);
end

% Optional title
title('Semicircular Polar Box Charts');
