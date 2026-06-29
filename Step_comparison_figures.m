clc; clear all; close all;

%% ========================================================================
%  1. Define path and set parameters
%  =======================================================================
workfolder = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\MatLab to use';
addpath(workfolder);

% Define the friction batch of interest
batches = {'Final_friction_45'};
chezy_values =[45];

% Create a main directory to save the figures
saveDirMain = fullfile(workfolder, 'Main_Text_Figures');
if ~exist(saveDirMain, 'dir')
    mkdir(saveDirMain);
end

% Physical Constants & Matrix parameters
g      = 9.81;            
Ts2    = 12 * 3600;       % S2 Period in seconds
w      = 2 * pi / Ts2;    % S2 Angular frequency [rad/s]
Qall   = 100;          
h0     = 1.00;         
lcAll  = 150000;          % 150 km mark where transition begins
Slopes = [0.1, 0.05, 0.02275, 0.02, 0.01137, 0.01, 0.0075, 0.005, 0.002]; 
dx     = 100;          

set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultLineLineWidth', 1.8);

% Preallocate error storage for Figure 3.4 
mean_error_seaward = zeros(length(batches), length(Slopes));

% Structure to hold results from all batches for plotting
AllBatchData = struct();

%% ========================================================================
%  2. MULTI-BATCH & MULTI-RUN PARAMETER LOOP
%  ========================================================================
for b = 1:length(batches)
    batch = batches{b};
    C_chezy = chezy_values(b);
    Cd = g / (C_chezy^2);
    
    baseDir = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\' batch filesep];
    fprintf('\n==================================================\n');
    fprintf('Processing Batch: %s (Chézy = %d)\n', batch, C_chezy);
    fprintf('==================================================\n');
    
    for pp = 1:length(Slopes)
        runName = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', abs(-Qall), h0*100, lcAll, pp);
        matFile = [baseDir runName filesep 'md.mat'];
        
        if ~isfile(matFile)
            warning('Data file not found: %s. Skipping...', matFile);
            continue;
        end
        
        load(matFile, 'md');
        Nx = size(md.x, 1);
        t_sec = md.time(1:end) * 60; 
        
        if ismatrix(md.x)
            x_m  = md.x(:, 2); zb_m = md.zb(:, 2);
        else
            x_m  = md.x;       zb_m = md.zb;
        end
        x_km = x_m / 1000;
        
        h_local = zb_m; 
        h_local(h_local < 0.1) = 0.1; 
        
        %% 3. Delft3D Extraction & Harmonic Analysis
        S2amp = zeros(Nx, 1); S4amp = zeros(Nx, 1);
        S2phase = zeros(Nx, 1); Ut = zeros(Nx, 1);
        
        for xx = 1:Nx
            eta0 = squeeze(md.wlv(:, xx, 2));
            [amp1, amp2, phase1, ~] = get_S2_S4(eta0, t_sec, w);
            S2amp(xx)   = amp1;   S4amp(xx)   = amp2;
            S2phase(xx) = phase1;
            
            if isfield(md, 'u')
                u0 = squeeze(md.v(:, xx, 2));
                [u_amp1, ~, ~, ~] = get_S2_S4(u0, t_sec, w);
                Ut(xx) = u_amp1;
            else
                Ut(xx) = 0.1; 
            end
        end
        
        phase2unwrapped = unwrap(S2phase);
        
        AllBatchData.(batch)(pp).x_km         = x_km;
        AllBatchData.(batch)(pp).h_local      = h_local;
        AllBatchData.(batch)(pp).Delft3D.kimag = gradient(phase2unwrapped, dx);
        AllBatchData.(batch)(pp).Delft3D.kreal = -gradient(log(S2amp), dx);
        
        %% 4. k0 Theory
        k0 = w ./ sqrt(g .* h_local);
        AllBatchData.(batch)(pp).k0.kimag = k0;
        AllBatchData.(batch)(pp).k0.kreal = zeros(Nx, 1);
        
        %% 5. kf Theory
        r = (8 / (3 * pi)) .* Cd .* Ut ./ h_local; 
        kf_complex = (w ./ sqrt(g .* h_local)) .* sqrt(1 - 1i .* r ./ w); 
        AllBatchData.(batch)(pp).kf.kimag = real(kf_complex);
        AllBatchData.(batch)(pp).kf.kreal = -imag(kf_complex);
        
        %% 6. kj Theory & Misfit
        dh_dx_analytical = zeros(Nx, 1);
        idx_start = find(x_m >= lcAll, 1, 'first');
        idx_shelf_ref = find(x_m >= 500000, 1, 'first');
        if isempty(idx_shelf_ref); shelf_depth = 6.0; else, shelf_depth = h_local(idx_shelf_ref); end
        
        idx_end = find(x_m >= lcAll & abs(h_local - shelf_depth) < 1e-3, 1, 'first');
        if isempty(idx_end); idx_end = find(x_m >= 700000, 1, 'first'); end
        
        delta_h = h_local(idx_end) - h_local(idx_start);
        delta_x = x_m(idx_end) - x_m(idx_start);
        if delta_x > 0
    analytical_slope = delta_h / delta_x;
else
    analytical_slope = 0;
end

        
        dh_dx_analytical(idx_start:idx_end) = analytical_slope;
        inv_La = -dh_dx_analytical ./ h_local; 
        if inv_La ~= 0
    La = 1 ./ inv_La;
else
    La = Inf;
end

        
        Delta = zeros(Nx, 1);
        valid_La = ~isinf(La) & (La ~= 0);
        Delta(valid_La) = sqrt(g .* h_local(valid_La)) ./ (2 .* La(valid_La) .* w); 
        
        kj_complex = (w ./ sqrt(g .* h_local)) .* sqrt(1 - Delta.^2 - 1i .* r ./ w); 
        kjimag = real(kj_complex);
        kjreal = -imag(kj_complex);
        
        AllBatchData.(batch)(pp).kj.kimag = kjimag;
        AllBatchData.(batch)(pp).kj.kreal = kjreal;
        
          % Calculate Relative Error
        k_j_mag   = sqrt(kjreal.^2 + kjimag.^2);
        k_mod_mag = sqrt(AllBatchData.(batch)(pp).Delft3D.kreal.^2 + AllBatchData.(batch)(pp).Delft3D.kimag.^2);
        epsilon_k = (abs(k_mod_mag - k_j_mag) ./ max(k_j_mag, 1e-6)) * 100;
        
        % Store mean seaward error (0 to 150 km zone)
        idx_seaward = x_m <= lcAll;
        mean_error_seaward(b, pp) = mean(epsilon_k(idx_seaward));
    end
end
%% ========================================================================
%  4. SECTION 3.5: STEP-LIKE REPRESENTATION OF LINEAR BATHYMETRIES
%  ========================================================================
fprintf('\n==================================================\n');
fprintf('Processing Discretization Experiments \n');
fprintf('==================================================\n');

% Define the stepped directories and their step counts
step_batches = {'Test_stairs_2', 'Test_stairs_4', 'Test_stairs_8'};
step_counts  = [2,4,8];

% Run 1 (Steepest/Sharp), Run 5 (Intermediate), Run 9 (Gentle Baseline)
slopes_to_compare =[1,2,4,5,8,9]; 

step_colors = {[0.85 0.33 0.1], [0.0 0.45 0.74], [0.47 0.67 0.19]}; 

for sl = 1:length(slopes_to_compare)
    pp = slopes_to_compare(sl);
    
    smooth_run = AllBatchData.Final_friction_45(pp);
    x_km = smooth_run.x_km;
    Nx = length(x_km);
   % Find where the transition zone begins (150 km)
idx_trans_start = find(x_km >= 150, 1, 'first');

% Look for where the depth matches 6.0 meters after the transition has started
% Look for where the absolute difference is very small (near zero)
idx_shelf_limit = find(x_km >= 150 & abs(smooth_run.h_local - 6.0) < 0.01, 1, 'first');

if ~isempty(idx_shelf_limit)
    % Set limit to where it hits the shelf plus a 30 km buffer
    x_limit_max = x_km(idx_shelf_limit) + 30;
else
    % Fallback if a run stabilizes at a different shelf depth (like 5.5m or 6.5m)
    % Find where the slope becomes completely flat again after 150 km
    idx_flat = find(x_km >= 150 & abs(gradient(smooth_run.h_local, dx)) < 1e-6, 1, 'first');
    if ~isempty(idx_flat)
        x_limit_max = x_km(idx_flat) + 30;
    else
        x_limit_max = 500; % Ultimate fallback layout limit
    end
end

% Ensure the maximum limit never tries to look past the end of the data domain
x_limit_max = min(x_limit_max, max(x_km));


    % Initialize a comparison structure for the stepped runs
    StepData = struct();
    
    % --- Step 1: Extract Delft3D Data from Step Folders ---
    for s = 1:length(step_batches)
        s_batch = step_batches{s};
        baseStepDir = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\' s_batch filesep];
        
        % Match run names consistently 
        runName = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', abs(-Qall), h0*100, lcAll, pp);
        matFile = [baseStepDir runName filesep 'md.mat'];
        
        if ~isfile(matFile)
            warning('Stepped run file not found: %s. Skipping...', matFile);
            continue;
        end
        
        % Load stepped output data structure
        step_load = load(matFile, 'md');
        
        % Process spatial parameters identically to section 2 loop
        if ismatrix(step_load.md.x)
            h_step_m = step_load.md.zb(:, 2);
        else
            h_step_m = step_load.md.zb;
        end
        h_step_m(h_step_m < 0.1) = 0.1;
        
        % Perform identical extraction and harmonic phase analysis loops
        S2amp_s = zeros(Nx, 1); S2phase_s = zeros(Nx, 1);
        t_sec_s = step_load.md.time(1:end) * 60;
        
        for xx = 1:Nx
            eta_s = squeeze(step_load.md.wlv(:, xx, 2));
            [amp1_s, ~, phase1_s, ~] = get_S2_S4(eta_s, t_sec_s, w);
            S2amp_s(xx)   = amp1_s;
            S2phase_s(xx) = phase1_s;
        end
        
        phase2unwrapped_s = unwrap(S2phase_s);
        
        % Detect step transitions by finding where the bathymetry changes sharply
is_boundary = abs(gradient(h_step_m, dx)) > 0.0005;

% Set the kreal and kimag values at these boundaries to NaN
kimag_clean = gradient(phase2unwrapped_s, dx);
kreal_clean = -gradient(log(S2amp_s), dx);

kimag_clean(is_boundary) = NaN;
kreal_clean(is_boundary) = NaN;

% Store the cleaned variables instead of the raw ones
StepData(s).kimag = kimag_clean;
StepData(s).kreal = kreal_clean;

        % Store derived properties for comparison plotting
        StepData(s).h_local = h_step_m;
        StepData(s).amp     = S2amp_s;
        StepData(s).kimag   = gradient(phase2unwrapped_s, dx);
        StepData(s).kreal   = -gradient(log(S2amp_s), dx);
    end
    
    % If data is missing for the run, do not throw errors
    if isempty(fieldnames(StepData)); continue; end
    
    % --- Step 2: Plot Comprehensive Hydrodynamic Panels (Amp, Phase, kimag, kreal) ---
    fig_panels = figure('Name', sprintf('Slope_%d_Step_Comparison', pp), ...
                    'Units', 'centimeters', ...
                    'Position', [2, 2, 16, 22]);
                set(gcf, 'Color', 'w');


    
    % Panel A: Bathymetry Profiles
    subplot(4, 1, 1); hold on; box on; grid on;
    plot(x_km, smooth_run.h_local, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Smooth Slope');
    for s = 1:length(step_counts)
        plot(x_km, StepData(s).h_local, 'Color', step_colors{s}, 'LineWidth', 1.5, 'DisplayName', sprintf('%d Steps', step_counts(s)));
    end
    xlim([0 x_limit_max]);
 ylim([0 14]); set(gca, 'YDir', 'reverse');
    ylabel('Depth (m)'); title(sprintf('Slope %d: Linear and Step-like Bathymetric Configurations', pp));
    legend('Location', 'NorthWest');
 slope_m_km_values = [1.00, 0.50, 0.20, 0.1137, 0.05, 0.02];
current_slope_val = slope_m_km_values(sl);

% Small slope annotation in top-right corner 
text(0.98, 0.95, sprintf('Slope %d: %.2f m/km', pp, current_slope_val), ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 8, ...
    'FontAngle', 'italic', ...
    'Color', [0.3 0.3 0.3], ...
    'BackgroundColor', 'w', ...
    'Margin', 1);


    
    % Panel B: Phase Propagation (kimag)
    subplot(4, 1, 2); hold on; box on; grid on;
    plot(x_km, smooth_run.Delft3D.kimag * 10^5, 'k-', 'LineWidth', 1.5);
    for s = 1:length(step_counts)
        plot(x_km, StepData(s).kimag * 10^5, 'Color', step_colors{s}, 'LineWidth', 1.5);

    end
    xlim([0 x_limit_max]);
 ylabel('k_{imag} (\times 10^{-5} rad/m)');
    title('Evolution of k_{imag}: Phase Propagation');
    
    % Panel C: Damping Exponent Energy Barriers (kreal)
    subplot(4, 1, 3); hold on; box on; grid on;
    plot(x_km, smooth_run.Delft3D.kreal * 10^5, 'k-', 'LineWidth', 1.5);
    for s = 1:length(step_counts)
        plot(x_km, StepData(s).kreal * 10^5, 'Color', step_colors{s}, 'LineWidth', 1.5);
    end
    xlim([0 x_limit_max]);
 ylabel('k_{real} (\times 10^{-5} rad/m)');
    title('Evolution of k_{real}: Amplitude Damping');
    
    % Panel D: Amplitude Modulation Fields
    subplot(4, 1, 4); hold on; box on; grid on;
    % Continuous profile data extraction for baseline
    t_sec = md.time(1:end) * 60;
    amp_smooth = zeros(Nx, 1);
t_sec = md.time(1:end) * 60;
for xx = 1:Nx
    eta0 = squeeze(md.wlv(:, xx, 2));
    [amp1, ~, ~, ~] = get_S2_S4(eta0, t_sec, w);
    amp_smooth(xx) = amp1;
end


    % Plot relative amplitude response comparisons
    for s = 1:length(step_counts)
        plot(x_km, StepData(s).amp, 'Color', step_colors{s}, 'LineWidth', 1.5);
    end
    xlim([0 x_limit_max]);
 xlabel('Distance along estuary (km)'); ylabel('Amplitude (m)');
    title('S_2 Amplitude Evolution');
    
    % Save consolidated manuscript matrix figures
    saveas(fig_panels, fullfile(saveDirMain, sprintf('Slope%d_Discretization_Panels.png', pp)));
end


%% ========================================================================
%  5. MISFIT ERROR SUMMARY PLOT 
%  ========================================================================
fprintf('\nGenerating Global Misfit Error Summary Figure...\n');

% Define the Delta values corresponding to your 6 specific evaluated slopes
% Mapping matches Table 1: Slopes = [0.1, 0.05, 0.02, 0.01137, 0.005, 0.002]
delta_values = [4.3957, 2.1979, 0.8791, 0.5, 0.2198, 0.0879];

% Preallocate error storage matrix: 3 step configurations x 6 evaluated slopes
misfit_summary_pct = zeros(length(step_batches), length(slopes_to_compare));

for sl = 1:length(slopes_to_compare)
    pp = slopes_to_compare(sl);
    smooth_run = AllBatchData.Final_friction_45(pp);
    
    idx_shelf_limit = find(smooth_run.x_km >= 150 & abs(smooth_run.h_local - 6.0) < 0.01, 1, 'first');
    if ~isempty(idx_shelf_limit)
        x_limit_max = smooth_run.x_km(idx_shelf_limit) + 30;
    else
        x_limit_max = 500;
    end
    x_limit_max = min(x_limit_max, max(smooth_run.x_km));
    active_idx = smooth_run.x_km <= x_limit_max;
    
    runName_control = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', abs(-Qall), h0*100, lcAll, pp);
    matFile_control = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\Final_friction_45' filesep runName_control filesep 'md.mat'];
    load(matFile_control, 'md');
    
    Nx_s = size(md.x, 1);
    amp_smooth_calc = zeros(Nx_s, 1);
    t_sec_c = md.time(1:end) * 60;
    for xx = 1:Nx_s
        eta_c = squeeze(md.wlv(:, xx, 2));
        [amp1_c, ~, ~, ~] = get_S2_S4(eta_c, t_sec_c, w);
        amp_smooth_calc(xx) = amp1_c;
    end
    
    for s = 1:length(step_batches)
        s_batch = step_batches{s};
        baseStepDir = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\' s_batch filesep];
        runName = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', abs(-Qall), h0*100, lcAll, pp);
        matFile = [baseStepDir runName filesep 'md.mat'];
        
        if isfile(matFile)
            step_load = load(matFile, 'md');
            S2amp_s = zeros(Nx_s, 1);
            t_sec_s = step_load.md.time(1:end) * 60;
            
            for xx = 1:Nx_s
                eta_s = squeeze(step_load.md.wlv(:, xx, 2));
                [amp1_s, ~, ~, ~] = get_S2_S4(eta_s, t_sec_s, w);
                S2amp_s(xx) = amp1_s;
            end
            
            % Calculate Mean Absolute Percentage Error inside active zone
            pct_diff = (abs(S2amp_s(active_idx) - amp_smooth_calc(active_idx)) ./ amp_smooth_calc(active_idx)) * 100;
            misfit_summary_pct(s, sl) = mean(pct_diff);
        end
    end
end

% Extract global maximum misfit values
max_err_2step = max(misfit_summary_pct(1, :));
max_err_4step = max(misfit_summary_pct(2, :));
max_err_8step = max(misfit_summary_pct(3, :));

% Create compact landscape figure
fig_misfit = figure('Name', 'Global_Step_Misfit_Percentage', 'Units', 'centimeters', 'Position',[2, 2, 16, 12]);
set(gcf, 'Color', 'w');
hold on; box on; grid on;

% Plot lines with basic legendary labels
plot(delta_values, misfit_summary_pct(1,:), 'o-', 'Color', step_colors{1}, 'LineWidth', 1.5, 'MarkerFaceColor', step_colors{1}, 'MarkerSize', 5, 'DisplayName', '2-Step Configuration');
plot(delta_values, misfit_summary_pct(2,:), 's-', 'Color', step_colors{2}, 'LineWidth', 1.5, 'MarkerFaceColor', step_colors{2}, 'MarkerSize', 5, 'DisplayName', '4-Step Configuration');
plot(delta_values, misfit_summary_pct(3,:), 'd-', 'Color', step_colors{3}, 'LineWidth', 1.5, 'MarkerFaceColor', step_colors{3}, 'MarkerSize', 5, 'DisplayName', '8-Step Configuration');

xlabel('Non-Dimensional Slope Parameter \Delta', 'FontSize', 10, 'FontName', 'Helvetica');
ylabel('Mean Amplitude Misfit Error (%)', 'FontSize', 10, 'FontName', 'Helvetica');
title('Discretization Percentage Error Scaling vs. Bed Geometry', 'FontSize', 11, 'FontName', 'Helvetica');


% Compile multiline string box for maximum values
max_error_string = sprintf(['Maximum Misfit Errors:\n' ...
                            '2-Step Representation: %.1f%%\n' ...
                            '4-Step Representation: %.1f%%\n' ...
                            '8-Step Representation: %.1f%%'], ...
                            max_err_2step, max_err_4step, max_err_8step);

% Keep the legend at the top-right corner
legend('Location', 'NorthEast');

% Place the text block directly underneath the legend box 
text(0.98, 0.65, max_error_string, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 8, ...
    'FontAngle', 'italic', ...
    'Color', [0.3 0.3 0.3], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);


saveas(fig_misfit, fullfile(saveDirMain, 'Global_Discretization_Misfit_Percentage.png'));
fprintf('Misfit percentage graph with text block annotation successfully saved.\n');


%% ========================================================================
%  6. SPATIAL MAXIMUM MISFIT ERROR PROFILES
%  =======================================================================
fprintf('\nGenerating Spatial Maximum Misfit Error Figure...\n');

% Set up a compact 2-panel figure (Top: Slope 1, Bottom: Slope 9)
fig_max_spatial = figure('Name', 'Spatial_Max_Misfit', 'Units', 'centimeters', 'Position', [2, 2, 16, 14]);
set(gcf, 'Color', 'w');

% Target the two extreme cases for clear comparison
target_slopes =[9,1]
plot_titles = {'Slope 1: Sharp Discretization Max Error Profile', 'Slope 9: Gentle Discretization Max Error Profile'};

for idx = 1:2
    pp = target_slopes(idx);
    smooth_run = AllBatchData.Final_friction_45(pp);
    
      idx_shelf_limit = find(smooth_run.x_km >= 150 & abs(smooth_run.h_local - 6.0) < 0.01, 1, 'first');
    if ~isempty(idx_shelf_limit)
        x_limit_max = smooth_run.x_km(idx_shelf_limit) + 30;
    else
        x_limit_max = 500;
    end
    x_limit_max = min(x_limit_max, max(smooth_run.x_km));
    active_idx = smooth_run.x_km <= x_limit_max;
    x_active = smooth_run.x_km(active_idx);
    
    % Load smooth control baseline amplitude
    runName_control = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', abs(-Qall), h0*100, lcAll, pp);
    matFile_control = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\Final_friction_45' filesep runName_control filesep 'md.mat'];
    load(matFile_control, 'md');
    
    Nx_s = size(md.x, 1);
    amp_smooth_calc = zeros(Nx_s, 1);
    t_sec_c = md.time(1:end) * 60;
    for xx = 1:Nx_s
        eta_c = squeeze(md.wlv(:, xx, 2));
        [amp1_c, ~, ~, ~] = get_S2_S4(eta_c, t_sec_c, w);
        amp_smooth_calc(xx) = amp1_c;
    end
    
    subplot(2, 1, idx); hold on; box on; grid on;
    
    % Extract and plot the spatial error profile for each step batch
    for s = 1:length(step_batches)
        s_batch = step_batches{s};
        baseStepDir = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\' s_batch filesep];
        runName = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', abs(-Qall), h0*100, lcAll, pp);
        matFile = [baseStepDir runName filesep 'md.mat'];
        
        if isfile(matFile)
            step_load = load(matFile, 'md');
            S2amp_s = zeros(Nx_s, 1);
            t_sec_s = step_load.md.time(1:end) * 60;
            
            for xx = 1:Nx_s
                eta_s = squeeze(step_load.md.wlv(:, xx, 2));
                [amp1_s, ~, ~, ~] = get_S2_S4(eta_s, t_sec_s, w);
                S2amp_s(xx) = amp1_s;
            end
            
            % Calculate local percentage error along every point of the channel
            local_pct_error = (abs(S2amp_s(active_idx) - amp_smooth_calc(active_idx)) ./ amp_smooth_calc(active_idx)) * 100;
            
            % Clean up boundary spikes 
            if isfield(StepData, 'h_local') || exist('h_step_m', 'var')
                % Use a loose threshold to keep the plot looking smooth but true
                is_boundary = abs(gradient(h_step_m, dx)) > 0.005; 
                local_pct_error(is_boundary(active_idx)) = NaN;
            end
            
            plot(x_active, local_pct_error, 'LineStyle', ':', 'Color', step_colors{s}, 'LineWidth', 1.2, 'DisplayName', sprintf('%d Steps', step_counts(s)));
        end
    end
    
    xlim([100 x_limit_max]); % Zoom explicitly into the transition region
    ylabel('Local Error (%)', 'FontSize', 9);
    title(plot_titles{idx}, 'FontSize', 10);
    if idx == 1; legend('Location', 'NorthWest'); end
    if idx == 2; xlabel('Distance along estuary (km)', 'FontSize', 10); end
end

saveas(fig_max_spatial, fullfile(saveDirMain, 'Spatial_Maximum_Misfit_Errors.png'));
fprintf('Spatial maximum error profiles saved successfully.\n');

