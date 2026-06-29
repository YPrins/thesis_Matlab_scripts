clc; clear; close all;

%% ============================================================
%  SETTINGS
%% ============================================================
workfolder = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\MatLab to use';
addpath(workfolder);

savefolder = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\Bathymetry_figure';

ls = 1200000;
bc = 500;
dx = 100;
dy = 500;

hd = [12, 6, -2];
lc = [150000, 700000];
Slopes = [0.1,0.05,0.02275,0.02,0.01137,0.01,0.0075,0.005,0.002];
colors = lines(length(Slopes));

%% ============================================================
%  FIGURE 
%% ============================================================
fig = figure('Color','w');
set(fig,'Units','pixels','Position',[100 100 1400 800]);

axMain = axes('Parent',fig);
hold(axMain,'on'); grid(axMain,'on');

%% ============================================================
%  MAIN PLOT
%% ============================================================
for pp = 1:length(Slopes)

    [depth, cor, ~, ~] = make_depth_trapjes( ...
        ls, dx, dy, bc, hd, lc, Slopes(pp));

    x_km = cor.x(1,:) / 1000;
    zb   = depth(1,:);

    lbl = sprintf('Run %d: S = %.3f (%.2f m/km)', ...
    pp, Slopes(pp), Slopes(pp)*10);
    

    plot(axMain, x_km, zb, ...
        'LineWidth', 2, ...
        'Color', colors(pp,:), ...
        'DisplayName', lbl);

end

set(axMain,'YDir','reverse');

xlabel(axMain,'Distance along estuary (km)');
ylabel(axMain,'Bed level (m)');
title(axMain,'Bathymetry with zoom on step region');

legendHandle = legend(axMain,'Location','southeast');
set(legendHandle,'Box','on');

%% ============================================================
%  STEP MARKER 
%% ============================================================
hStep = xline(axMain,150,'--k','Step location', ...
    'LabelVerticalAlignment','top', ...
    'LineWidth',1.2);

hStep.Annotation.LegendInformation.IconDisplayStyle = 'off';
%% ============================================================
%  FORCE MAIN AXES TO BACK
%% ============================================================
axMain.Layer = 'bottom';

%% ============================================================
%  INSET 
%% ============================================================
axInset = axes('Parent',fig, ...
    'Position',[0.58 0.58 0.32 0.32], ...
    'Box','on');

hold(axInset,'on'); grid(axInset,'on');

for pp = 1:length(Slopes)

    [depth, cor, ~, ~] = make_depth_trapjes( ...
        ls, dx, dy, bc, hd, lc, Slopes(pp));

    x_km = cor.x(1,:) / 1000;
    zb   = depth(1,:);

    idx = x_km >= 140 & x_km <= 300;

    plot(axInset, x_km(idx), zb(idx), ...
        'LineWidth', 1.5, ...
        'Color', colors(pp,:));
end

set(axInset,'YDir','reverse');
xlim(axInset,[140 300]);
ylim(axInset,[5.5 12.5]);

title(axInset,'Zoom: Step region');
xlabel(axInset,'km');
ylabel(axInset,'m');

axInset.Layer = 'top';

drawnow;  

uistack(axMain,'bottom');
uistack(axInset,'top');
uistack(legendHandle,'top');
set(findall(fig,'Type','line'),'HandleVisibility','on');

if ~exist(savefolder,'dir')
    mkdir(savefolder);
end

% PNG high resolution
print(fig, fullfile(savefolder,'bathymetry_with_inset.png'), '-dpng', '-r300');