clc; clear all; close all;

%% ========================================================================
%  1. Define path and set parameters
%  ========================================================================
workfolder = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\MatLab to use';
batch      = 'Final_friction_45';  % Change this and the friction value together! Change 3 times!!
baseDir    = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\' batch filesep];

addpath(workfolder);

% Create automatic save directory for figures based on friction parameter
saveDir = fullfile(workfolder, sprintf('Visualisation friction %d', 45)); % Update with C_chezy matching value
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Physical Constants
g = 9.81;            
C_chezy = 45;          % Change this value !!!
Cd = g / (C_chezy^2);  % Drag coefficient definition

% Tidal Settings (S2 forcing)
Ts2 = 12 * 3600;       % S2 Period in seconds (12 hours exactly)
w   = 2 * pi / Ts2;    % S2 Angular frequency [rad/s]

% Matrix parameters
Qall   = 100;          % Discharge (m3/s)
h0     = 1.00;         % Tidal amplitude (m)
lcAll  = 150000;       % Characteristic length (150 km (where the step starts))
Slopes = [0.1,0.05,0.02275,0.02,0.01137,0.01,0.0075,0.005,0.002]; % bed slopes
%[0.1,0.05,0.02275,0.02,0.01137,0.01,0.0075,0.005,0.002]; For linear runs
%[0.1,0.05,0.02,0.01137,0.005,0.002]; For stair runs

dx     = 100;          % Grid resolution in meters

% Global Graphics Theme Configuration
set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultLineLineWidth', 1.8);

% Initialize centralized figures 
figPlot1 = figure('Name', 'Plot 1: Global Hydrodynamic Response', 'Color', 'w', 'Position', [50, 50, 950, 950]);
figPlot3 = figure('Name', 'Plot 3: Analytical Misfit Relative Error', 'Color', 'w', 'Position', [150, 50, 950, 700]);

% Custom Thesis Color Scheme Matrix
colors_matrix = [
    0.0000 0.4470 0.7410;  % blue
    0.8500 0.3250 0.0980;  % orange
    0.9290 0.6940 0.1250;  % yellow
    0.4940 0.1840 0.5560;  % purple
    0.4660 0.6740 0.1880;  % green
    0.3010 0.7450 0.9330;  % cyan
    0.6350 0.0780 0.1840;  % dark red
    0.0000 0.0000 0.0000;  % black
    0.7500 0.7500 0.7500;  % grey
];

Results = struct();

h_leg_p1 = zeros(length(Slopes), 1);
h_leg_p3 = zeros(length(Slopes), 1);

%% ========================================================================
%  2. MULTI-RUN PARAMETER LOOP
%  ========================================================================
for pp = 1:length(Slopes)
    
    current_color = colors_matrix(pp, :);
    
    runName = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', abs(-Qall), h0*100, lcAll, pp);
    matFile = [baseDir runName filesep 'md.mat'];
    
    if ~isfile(matFile)
        warning('Data file not found for processing: %s. Skipping...', matFile);
        continue;
    end
    
    load(matFile, 'md');
    Nx = size(md.x, 1);
    t_sec = md.time(1:end) * 60; 
    
    if ismatrix(md.x)
        x_m  = md.x(:, 2); 
        zb_m = md.zb(:, 2);
    else
        x_m  = md.x;
        zb_m = md.zb;
    end
    x_km = x_m / 1000;
    
    % True positive water depth profile h(x)
    h_local = zb_m; 
    h_local(h_local < 0.1) = 0.1; % Floor boundary to prevent division by zero
    
    %% --------------------------------------------------------------------
    %  3. Delft3D Extraction & Harmonic Analysis 
    %  --------------------------------------------------------------------
    S2amp   = zeros(Nx, 1); S4amp   = zeros(Nx, 1);
    S2phase = zeros(Nx, 1); S4phase = zeros(Nx, 1);
    Ut      = zeros(Nx, 1);
    
    for xx = 1:Nx
        eta0 = squeeze(md.wlv(:, xx, 2));
        [amp1, amp2, phase1, phase2] = get_S2_S4(eta0, t_sec, w);
        
        S2amp(xx)   = amp1;   S4amp(xx)   = amp2;
        S2phase(xx) = phase1; S4phase(xx) = phase2;
        
        if isfield(md, 'u')
            u0 = squeeze(md.v(:, xx, 2));
            [u_amp1, ~, ~, ~] = get_S2_S4(u0, t_sec, w);
            Ut(xx) = u_amp1;
        else
            Ut(xx) = 0.1; 
        end
    end
    
    phase2unwrapped = unwrap(S2phase);
    
    Results(pp).Delft3D.S2amp    = S2amp;
    Results(pp).Delft3D.S4amp    = S4amp;
    Results(pp).Delft3D.S2phase  = phase2unwrapped;
    Results(pp).Delft3D.S4phase  = unwrap(S4phase);
    Results(pp).Delft3D.kimag    = gradient(phase2unwrapped, dx);
    Results(pp).Delft3D.kreal    = -gradient(log(S2amp), dx);
    Results(pp).Delft3D.celerity = w ./ max(Results(pp).Delft3D.kimag, 1e-6);
    
    %% --------------------------------------------------------------------
    %  4. k0 Theory - Integrated Phase
    %  --------------------------------------------------------------------
    c0 = sqrt(g .* h_local);
    k0 = w ./ sqrt(g .* h_local);
    
    Results(pp).k0.kreal     = zeros(Nx, 1);
    Results(pp).k0.kimag     = k0;
    Results(pp).k0.celerity  = c0;
    Results(pp).k0.phase     = cumtrapz(x_m, k0); 
    Results(pp).k0.amplitude = S2amp(1) .* ones(Nx, 1); 
    
    %% --------------------------------------------------------------------
    %  5. kf Theory - Integrated Phase & Damping
    %  --------------------------------------------------------------------
    r = (8 / (3 * pi)) .* Cd .* Ut ./ h_local; 
    kf_complex = (w ./ sqrt(g .* h_local)) .* sqrt(1 - 1i .* r ./ w); 
    
    kfimag = real(kf_complex);
    kfreal = -imag(kf_complex);
    
    Results(pp).kf.kreal     = kfreal;
    Results(pp).kf.kimag     = kfimag;
    Results(pp).kf.celerity  = w ./ max(kfimag, 1e-6);
    Results(pp).kf.phase     = cumtrapz(x_m, kfimag);
    Results(pp).kf.amplitude = S2amp(1) .* exp(-cumtrapz(x_m, kfreal));
    
    %% --------------------------------------------------------------------
    %  6. kj Theory & Relative Error 
    %  --------------------------------------------------------------------
    dh_dx_analytical = zeros(Nx, 1);
    
    idx_start = find(x_m >= lcAll, 1, 'first');
    idx_shelf_ref = find(x_m >= 500000, 1, 'first');
    if isempty(idx_shelf_ref)
        shelf_depth = 6.0; 
    else
        shelf_depth = h_local(idx_shelf_ref);
    end
    
    idx_end = find(x_m >= lcAll & abs(h_local - shelf_depth) < 1e-3, 1, 'first');
    if isempty(idx_end)
        idx_end = find(x_m >= 700000, 1, 'first'); 
    end
    
    idx_before_step = (1:Nx)' < idx_start;
    idx_step        = (1:Nx)' >= idx_start & (1:Nx)' <= idx_end;
    idx_after_step  = (1:Nx)' > idx_end;
    
    delta_h = h_local(idx_end) - h_local(idx_start);
    delta_x = x_m(idx_end) - x_m(idx_start);
    
    if delta_x > 0
        analytical_slope = delta_h / delta_x;
    else
        analytical_slope = 0;
    end
    
    dh_dx_analytical(idx_before_step) = 0;
    dh_dx_analytical(idx_step)        = analytical_slope;
    dh_dx_analytical(idx_after_step)  = 0;
    
    inv_La = -dh_dx_analytical ./ h_local; 
    La = zeros(Nx, 1);
    La(inv_La ~= 0) = 1 ./ inv_La(inv_La ~= 0);
    La(inv_La == 0) = Inf; 
    
    Delta = zeros(Nx, 1);
    valid_La = ~isinf(La) & (La ~= 0);
    Delta(valid_La) = sqrt(g .* h_local(valid_La)) ./ (2 .* La(valid_La) .* w); 
    Delta(~valid_La) = 0; 
    
    kj_complex = (w ./ sqrt(g .* h_local)) .* sqrt(1 - Delta.^2 - 1i .* r ./ w); 
    kjimag = real(kj_complex);
    kjreal = -imag(kj_complex);
    
    Results(pp).kj.kreal     = kjreal;
    Results(pp).kj.kimag     = kjimag;
    Results(pp).kj.celerity  = w ./ max(kjimag, 1e-6);
    Results(pp).kj.phase     = cumtrapz(x_m, kjimag);
    Results(pp).kj.amplitude = S2amp(1) .* exp(-cumtrapz(x_m, kjreal));
    
    k_j_mag   = sqrt(kjreal.^2 + kjimag.^2);
    k_mod_mag = sqrt(Results(pp).Delft3D.kreal.^2 + Results(pp).Delft3D.kimag.^2);
    Results(pp).RelativeError = (abs(k_mod_mag - k_j_mag) ./ max(k_j_mag, 1e-6)) * 100;

    %% --------------------------------------------------------------------
    %  6b. DEBUGGING DISPLAY
    %  --------------------------------------------------------------------
    max_Delta_in_step = max(Delta(idx_step));
    step_length_km = delta_x / 1000;
    
    mid_idx = round((idx_start + idx_end) / 2);
    if analytical_slope ~= 0
        Lh_km = h_local(mid_idx) / abs(analytical_slope) / 1000;
    else
        Lh_km = Inf;
    end
    
    fprintf('Run %2d | Slope Input: %7.5f | Step Length: %6.1f km | Lh: %5.1f km | Max Delta: %6.4f | 2Delta: %6.4f\n', ...
            pp, Slopes(pp), step_length_km, Lh_km, max_Delta_in_step, 2 * max_Delta_in_step);

    %% ========================================================================
    %  7. POPULATE PLOT 1 
    %  ========================================================================
    figure(figPlot1);
    label_text = sprintf('Slope %d: %.3f m/km', pp, Slopes(pp) * 10);
    
    % Panel C: S2 Amplitude (Row 1 of 3)
    subplot(3, 1, 1); hold on;
    h_leg_p1(pp) = plot(x_km, Results(pp).Delft3D.S2amp, 'Color', colors_matrix(pp,:), 'DisplayName', label_text);
    
    % Panel D: S4 Amplitude (Row 2 of 3)
    subplot(3, 1, 2); hold on;
    plot(x_km, Results(pp).Delft3D.S4amp, 'Color', colors_matrix(pp,:));
    
    % Panel E: S2 Spatial Phase (Row 3 of 3)
    subplot(3, 1, 3); hold on;
    plot(x_km, rad2deg(Results(pp).Delft3D.S2phase), 'Color', colors_matrix(pp,:));
    
    %% ========================================================================
    %  8. POPULATE PLOT 3 (Relative Error Master)
    %  ========================================================================
    figure(figPlot3);
    
    subplot(2, 1, 1); hold on;
    h_leg_p3(pp) = plot(x_km, h_local, 'Color', colors_matrix(pp,:), 'DisplayName', label_text);
    
    subplot(2, 1, 2); hold on;
    plot(x_km, Results(pp).RelativeError, 'Color', colors_matrix(pp,:));
    
    %% ========================================================================
    %  9. PLOT 2 (Per Run Breakdown)
    %  ========================================================================
    figDiag = figure('Name', sprintf('Plot 2: Run %d Diagnostics', pp), 'Color', 'w', 'Position', [200, 100, 850, 750]);
    
    subplot(3, 1, 1);
    plot(x_km, h_local, 'k-', 'LineWidth', 2);
    ylabel('Depth (m)'); title(sprintf('Run %d: Bathymetry Profile', pp), 'FontWeight', 'bold'); grid on;
    set(gca, 'YDir', 'reverse', 'Box', 'off'); xlim([0 500]); ylim([0 13]);
    
    subplot(3, 1, 2); hold on;
    plot(x_km, Results(pp).Delft3D.kimag, 'b-', 'LineWidth', 2, 'DisplayName', 'Delft3D');
    plot(x_km, Results(pp).k0.kimag, 'g--', 'DisplayName', 'k_0 (Frictionless)');
    plot(x_km, Results(pp).kf.kimag, 'r-.', 'DisplayName', 'k_f (Frictional)');
    plot(x_km, Results(pp).kj.kimag, 'm:', 'LineWidth', 2, 'DisplayName', 'k_j (Jay)');
    ylabel('k_{imag} (rad/m)'); title('Wavenumber Mechanics: Phase Propagation (k_{imag})', 'FontWeight', 'bold'); grid on;
    xlim([0 500]); set(gca, 'Box', 'off');
    kimag_step_max = max(Results(pp).kj.kimag(x_km <= 500));
    ylim([0, max(kimag_step_max * 2.5, 5e-5)]); 
    legend('show', 'Location', 'best', 'Box', 'off');
    
    subplot(3, 1, 3); hold on;
    plot(x_km, Results(pp).Delft3D.kreal, 'b-', 'LineWidth', 2);
    plot(x_km, Results(pp).k0.kreal, 'g--');
    plot(x_km, Results(pp).kf.kreal, 'r-.');
    plot(x_km, Results(pp).kj.kreal, 'm:', 'LineWidth', 2);
    ylabel('k_{real} (rad/m)'); xlabel('Distance along estuary (km)'); title('Wavenumber Mechanics: Damping Exponent (k_{real})', 'FontWeight', 'bold'); grid on;
    xlim([0 500]); set(gca, 'Box', 'off');
    
    idx_500 = x_km <= 500;
    max_d3d = max(abs(Results(pp).Delft3D.kreal(idx_500)));
    max_jay = max(abs(Results(pp).kj.kreal(idx_500)));
    kreal_visual_max = max(max_d3d, max_jay);
    ylim([-kreal_visual_max * 0.2, kreal_visual_max * 1.2]);
    
    saveas(figDiag, fullfile(saveDir, sprintf('Run_%d_Diagnostics.png', pp)));
    
    %% --------------------------------------------------------------------
    %  9b. Calculate Single Average Relative Error (0-150 km)
    %  --------------------------------------------------------------------
    idx_zone = (x_km >= 0) & (x_km <= 150);
    Results(pp).AverageRelativeError = mean(Results(pp).RelativeError(idx_zone));
    
    fprintf('Run %2d Summary | Average Reflection Misfit Error (0-150km): %5.2f%%\n', ...
            pp, Results(pp).AverageRelativeError);
end

%% ========================================================================
%  10. Polish Global Overlays & Save Master Charts 
%  ========================================================================

% Format Plot 1 (Global Hydrodynamics - 3 Subplots) 
figure(figPlot1);
subplot(3,1,1); ylabel('Amplitude (m)'); title('Panel C: S_2 Free-Surface Water Level Amplitude Profile', 'FontWeight', 'bold'); grid on; set(gca, 'Box', 'off'); xlim([0 500]); 
subplot(3,1,2); ylabel('Amplitude (m)'); title('Panel D: S_4 Higher Harmonic Amplitude Profile', 'FontWeight', 'bold'); grid on; set(gca, 'Box', 'off'); xlim([0 500]);
subplot(3,1,3); ylabel('Phase (deg)'); xlabel('Distance along estuary (km)'); title('Panel E: S_2 Spatial Phase Shift Profiles', 'FontWeight', 'bold'); grid on; set(gca, 'Box', 'off'); xlim([0 500]);

drawnow;
print(figPlot1, fullfile(saveDir, 'Global_Hydrodynamic_Response'), '-dpng', '-r300');

% Format Plot 3 (Spatial Misfit Profile)
figure(figPlot3);
subplot(2,1,1); ylabel('Depth (m)'); title('Panel A: Bathymetric Configurations', 'FontWeight', 'bold'); grid on; set(gca, 'YDir', 'reverse', 'Box', 'off'); xlim([0 500]); ylim([0 13]); legend(h_leg_p3, 'Location', 'eastoutside');
subplot(2,1,2); ylabel('Relative Error \epsilon_k (%)'); xlabel('Distance along estuary (km)'); title('Panel B: Analytical Misfit Error (\epsilon_k) Evolution', 'FontWeight', 'bold'); grid on; set(gca, 'Box', 'off'); 

xlim([0 300]); 
error_step_max = 0;
for pp = 1:length(Slopes)
    if isfield(Results(pp), 'RelativeError') && ~isempty(Results(pp).RelativeError)
        current_max = max(Results(pp).RelativeError(x_km <= 300));
        if current_max > error_step_max, error_step_max = current_max; end
    end
end
ylim([0, max(error_step_max * 1.5, 40)]); 

y_lims = ylim;
patch([0 150 150 0], [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], [0.9 0.9 0.9], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(75, y_lims(2)*0.85, 'Seaward Of Onset of Slope', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);

drawnow;
print(figPlot3, fullfile(saveDir, 'Analytical_Misfit_Relative_Error'), '-dpng', '-r300');

%% ========================================================================
%  11. Plot 4 - The Reflection Threshold Chart (Single Value Summary)
%  ========================================================================
avg_errors = zeros(length(Slopes), 1);
for pp = 1:length(Slopes)
    avg_errors(pp) = Results(pp).AverageRelativeError;
end

figPlot4 = figure('Name', 'Plot 4: Single Value Reflection Misfit Summary', 'Color', 'w', 'Position', [150, 150, 900, 600]);

hold on;
plot(Slopes, avg_errors, 'o-', 'Color', colors_matrix(1,:), 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', colors_matrix(1,:));

set(gca, 'XDir', 'reverse', 'XScale', 'log', 'Box', 'off'); 
xlabel('Bed Slope (m/km)');
ylabel('Mean Seaward Misfit Error \epsilon_k (%) [0-150 km]');
title(sprintf('Analytical Model Break-Down vs. Bed Steepness (Chézy = %d)', C_chezy), 'FontWeight', 'bold');
grid on;

text(0.02, 0.05, 'Steep slope', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
text(0.98, 0.05, 'Gentle slope', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'right');

for pp = 1:length(Slopes)
    text(Slopes(pp), avg_errors(pp) + max(avg_errors)*0.03, sprintf('Run %d', pp), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);
end

drawnow;
print(figPlot4, fullfile(saveDir, 'Reflection_Threshold_Misfit_Summary'), '-dpng', '-r300');

