clc; clear all; close all;

%% ========================================================================
%  1. CONFIGURATION & BATCH LIST
%  ========================================================================
workfolder = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\MatLab to use';
addpath(workfolder);

% Define all friction configurations to build the comparative figures
batches    = {'Final_friction_30', 'Final_friction_45', 'Final_friction_60'};
chezy_vals = [30, 45, 60];

Slopes = [0.1, 0.05, 0.02275, 0.02, 0.01137, 0.01, 0.0075, 0.005, 0.002]; 
g      = 9.81;            
Ts2    = 12 * 3600;       
w      = 2 * pi / Ts2;    
dx     = 100;
lcAll  = 150000;

% Initialize data storage for cross-batch comparison
SummaryError = zeros(length(Slopes), length(batches));
SpatialData  = struct();

%% ========================================================================
%  2. CORE PROCESSING LOOP (ALL BATCHES x ALL SLOPES)
%  ========================================================================
for b = 1:length(batches)
    batch   = batches{b};
    C_chezy = chezy_vals(b);
    Cd      = g / (C_chezy^2);
    baseDir = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\' batch filesep];
    
    fprintf('Processing Batch: %s (Chezy %d)...\n', batch, C_chezy);
    
    for pp = 1:length(Slopes)
        runName = sprintf('run_Qr%d_h0100_zb1%d_slope%d', 100, lcAll, pp);
        matFile = [baseDir runName filesep 'md.mat'];
        
        if ~isfile(matFile), continue; end
        load(matFile, 'md');
        
        Nx = size(md.x, 1);
        t_sec = md.time * 60;
        x_m  = md.x(:, min(2, size(md.x,2)));
        zb_m = md.zb(:, min(2, size(md.zb,2)));
        x_km = x_m / 1000;
        h_local = max(zb_m, 0.1);
        
        % --- 3. Extract Delft3D Data ---
        S2amp = zeros(Nx,1); S2phase = zeros(Nx,1); Ut = zeros(Nx,1);
        for xx = 1:Nx
            [amp1, ~, phase1, ~] = get_S2_S4(squeeze(md.wlv(:, xx, 2)), t_sec, w);
            S2amp(xx)   = amp1;
            S2phase(xx) = phase1;
            if isfield(md, 'u')
                [u_amp1, ~, ~, ~] = get_S2_S4(squeeze(md.v(:, xx, 2)), t_sec, w);
                Ut(xx) = u_amp1;
            else, Ut(xx) = 0.1; end
        end
        phase_unwrapped = unwrap(S2phase);
        kreal_mod = -gradient(log(S2amp), dx);
        kimag_mod = gradient(phase_unwrapped, dx);
        
        % --- 4. Analytical k_j Calculation ---
        r = (8 / (3 * pi)) .* Cd .* Ut ./ h_local;
        dh_dx = zeros(Nx, 1);
        idx_start = find(x_m >= lcAll, 1, 'first');
        idx_shelf = find(x_m >= 500000, 1, 'first');
        shelf_depth = h_local(idx_shelf);
        idx_end = find(x_m >= lcAll & abs(h_local - shelf_depth) < 1e-3, 1, 'first');
        
        if ~isempty(idx_start) && ~isempty(idx_end)
            dh_dx(idx_start:idx_end) = (h_local(idx_end) - h_local(idx_start)) / (x_m(idx_end) - x_m(idx_start));
        end
        
        inv_La = -dh_dx ./ h_local;
        Delta = zeros(Nx,1);
        idx_v = inv_La ~= 0;
        Delta(idx_v) = sqrt(g .* h_local(idx_v)) ./ (2 .* (1./inv_La(idx_v)) .* w);
        
        kj_complex = (w ./ sqrt(g .* h_local)) .* sqrt(1 - Delta.^2 - 1i .* r ./ w);
        kjreal = -imag(kj_complex);
        kjimag = real(kj_complex);
        
        % --- 5. Calculate Misfit Error (Seaward Reflection Zone 0-150km) ---
        k_j_mag   = sqrt(kjreal.^2 + kjimag.^2);
        k_mod_mag = sqrt(kreal_mod.^2 + kimag_mod.^2);
        rel_err   = (abs(k_mod_mag - k_j_mag) ./ max(k_j_mag, 1e-6)) * 100;
        
        % Averaging inside the flat seaward zone up to the foot of the step
        idx_zone = (x_km >= 0) & (x_km <= 150);
        SummaryError(pp, b) = mean(rel_err(idx_zone));
        
        % Save specific data for spatial profiling figures (using standard Chézy 45 as baseline)
        if C_chezy == 45
            SpatialData(pp).x_km    = x_km;
            SpatialData(pp).S2amp   = S2amp;
            SpatialData(pp).S2phase = phase_unwrapped;
            SpatialData(pp).kjreal  = kjreal;
            SpatialData(pp).kjimag  = kjimag;
            SpatialData(pp).kreal_m = kreal_mod;
            SpatialData(pp).kimag_m = kimag_mod;
        end
    end
end

%% ========================================================================
%  3. GRAPHICS GENERATION 
%  ========================================================================
set(0, 'DefaultAxesFontName', 'Helvetica', 'DefaultAxesFontSize', 11);
set(0, 'DefaultLineLineWidth', 1.5);

%% --- FIGURE 1: Global Hydrodynamic Spatial Response ---
fig1 = figure('Color', 'w', 'Position', [100, 100, 800, 600]); 
cmap = [0.8 0 0; 0 0.5 0; 0 0 0.8]; 
target_slopes =[1, 3, 9]; % Indices for Slope 1, Slope 3, and Slope 9
labels = {'Slope 1: Strong Funnell/Reflection', 'Slope 3: Critical Convergence', 'Slope 9: Friction Dominated'};

for i = 1:3
    pp = target_slopes(i);
    if ~isfield(SpatialData, 'x_km') || length(SpatialData) < pp || isempty(SpatialData(pp).x_km), continue; end
    
    % Subplot A: S2 Amplitude Profile
    subplot(2,1,1); hold on;
    plot(SpatialData(pp).x_km, SpatialData(pp).S2amp, 'Color', cmap(i,:), 'DisplayName', labels{i});
    title('Longitudinal Estuary Tidal Hydrodynamics (Chézy = 45)');
    ylabel('S_2 Amplitude (m)'); grid on; xlim([0 400]);
    legend('Location', 'northeast', 'Box', 'off');
    
    % Subplot B: Unwrapped Phase Profile
    subplot(2,1,2); hold on;
    plot(SpatialData(pp).x_km, SpatialData(pp).S2phase, 'Color', cmap(i,:));
    xlabel('Distance along estuary (km)'); ylabel('S_2 Phase Shift (rad)'); grid on; xlim([0 400]);
end

%% --- FIGURE 2: Analytical Misfit Breakdown (0 to 200 km) ---
fig2 = figure('Color', 'w', 'Position', [150, 150, 1000, 450]); 
target_comparison =[2, 7]; % Slope 2 (Abrupt) vs Slope 7 (Gentle)

for plot_idx = 1:2
    pp = target_comparison(plot_idx);
    if ~isfield(SpatialData, 'x_km') || length(SpatialData) < pp || isempty(SpatialData(pp).x_km), continue; end
    
    subplot(1, 2, plot_idx); hold on;
    plot(SpatialData(pp).x_km, SpatialData(pp).kimag_m, 'b-', 'DisplayName', 'Delft3D k_{imag}');
    plot(SpatialData(pp).x_km, SpatialData(pp).kjimag, 'b--', 'DisplayName', 'Jay Theory k_{imag}');
    plot(SpatialData(pp).x_km, SpatialData(pp).kreal_m, 'r-', 'DisplayName', 'Delft3D k_{real}');
    plot(SpatialData(pp).x_km, SpatialData(pp).kjreal, 'r--', 'DisplayName', 'Jay Theory k_{real}');
    
    grid on; xlim([0 200]); % Focused tightly on the seaward boundary + step entrance
    xlabel('Distance along estuary (km)'); ylabel('Wavenumber Magnitude (m^{-1})');
    title(sprintf('Wavenumber Mechanics: Run %d', pp));
    if plot_idx == 1, legend('Location', 'best', 'Box', 'off'); end
end

%% --- FIGURE 3: Friction vs Reflection Threshold  ---
fig3 = figure('Color', 'w', 'Position', [200, 200, 700, 500]);
hold on;
plot(Slopes, SummaryError(:,1), 'o-', 'Color', [0.85 0.33 0.1], 'MarkerFaceColor', [0.85 0.33 0.1]);
plot(Slopes, SummaryError(:,2), 's-', 'Color', [0.47 0.67 0.19], 'MarkerFaceColor', [0.47 0.67 0.19]);
plot(Slopes, SummaryError(:,3), 'd-', 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);

set(gca, 'XDir', 'reverse', 'XScale', 'log'); 
xlabel('Transition Bed Slope (m/km) \leftarrow [Gentle --- Steep]');
ylabel('Mean Seaward Misfit Error \epsilon_k (%) in 0-150 km Zone');
title('Reflection-Induced Error Growth Before the Bathymetric Step');
grid on;
legend({'Chézy 30 (High Friction)', 'Chézy 45 (Standard)', 'Chézy 60 (Low Friction)'}, 'Location', 'Northwest', 'Box', 'off');

y_limits = ylim;
patch([0.02275 0.1 0.1 0.02275], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
      [0.9 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.04, y_limits(2)*0.8, 'Reflection Breakout Zone', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');