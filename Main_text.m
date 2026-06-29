clc; clear all; close all;

%% ========================================================================
%  1. Define path and set parameters
%  =======================================================================
workfolder = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\MatLab to use';
addpath(workfolder);

% Define the three friction batches to loop through
batches = {'Final_friction_30', 'Final_friction_45', 'Final_friction_60'};
chezy_values =[30, 45, 60];

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

% Preallocate error storage for Figure 3.4 (3 friction cases x 9 slopes)
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
        analytical_slope = idxchoose(delta_x > 0, delta_h / delta_x, 0);
        
        dh_dx_analytical(idx_start:idx_end) = analytical_slope;
        inv_La = -dh_dx_analytical ./ h_local; 
        La = idxchoose(inv_La ~= 0, 1 ./ inv_La, Inf);
        
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
%  3. FIGURE GENERATION SECTION
%  ========================================================================
showcase_batch = 'Final_friction_45'; % Change for friction value

%% --- FIGURE 3.1: Progressive Breakdown Matrix ---
fig31 = figure('Name', 'Figure 3.1: Analytical Breakdown Matrix', 'Color', 'w', 'Position', [50, 50, 1500, 750]);
target_runs = [9, 5, 4, 1]; % Run 9 (Gentle), Run 5 (Intermediate), Run 4 (Sharp), Run 1 (Steep)
run_labels = {'Slope 9 (Gentle Slope)', 'Slope 5 (Intermediate Slope)', 'Slope 4 (Sharp Slope)', 'Slope 1 (Steep Step)'};
slope_texts = {
    'Slope 9: 0.02 m/km'
    'Slope 5: 0.1137 m/km'
    'Slope 4: 0.20 m/km'
    'Slope 1: 1.00 m/km'
};
for col = 1:4
    pp = target_runs(col);
    x_km = AllBatchData.(showcase_batch)(pp).x_km;
    
  % Row 1: kimag (Wave Propagation)
subplot(2, 4, col); hold on; grid on;

plot(x_km, AllBatchData.(showcase_batch)(pp).Delft3D.kimag * 1e5, 'b', 'LineWidth', 2);
plot(x_km, AllBatchData.(showcase_batch)(pp).k0.kimag * 1e5, 'g--');
plot(x_km, AllBatchData.(showcase_batch)(pp).kf.kimag * 1e5, 'r-.');
plot(x_km, AllBatchData.(showcase_batch)(pp).kj.kimag * 1e5, 'm:');

xlim([0, 300]);
ylim([0, 3.5]);

title(run_labels{col}, 'FontSize', 11, 'FontWeight', 'bold');

% Small slope annotation in lower-right corner
text(0.98, 0.05, slope_texts{col}, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 8, ...
    'FontAngle', 'italic', ...
    'Color', [0.3 0.3 0.3], ...
    'BackgroundColor', 'w', ...
    'Margin', 1);

if col == 1
    ylabel('k_{imag} (\times10^{-5} rad/m) (phase shift)', ...
        'FontWeight', 'bold');
end

if col == 1
    legend('Delft3D', ...
           'k_0 (Frictionless)', ...
           'k_f (Frictional)', ...
           'k_j (Jay)', ...
           'Location', 'northwest');
end
    
    % Row 2: kreal (Wave Damping)
    subplot(2, 4, col + 4); hold on; grid on; 
    plot(x_km, AllBatchData.(showcase_batch)(pp).Delft3D.kreal * 1e5, 'b', 'LineWidth', 2);
    plot(x_km, AllBatchData.(showcase_batch)(pp).k0.kreal * 1e5, 'g--');
    plot(x_km, AllBatchData.(showcase_batch)(pp).kf.kreal * 1e5, 'r-.');
    plot(x_km, AllBatchData.(showcase_batch)(pp).kj.kreal * 1e5, 'm:');
    xlim([0, 300]); 
    xlabel('Distance along estuary (km)', 'FontWeight', 'bold');
    
    if col == 1
        ylabel('k_{real} (\times10^{-5} rad/m) (amplitude damping)', 'FontWeight', 'bold'); 
    end
end
saveas(fig31, fullfile(saveDirMain, 'Figure_3_1_Breakdown_Matrix.png'));
%% --- FIGURE 3.2: Spatial Phase Evolution (kimag Profiles) ---
fig32 = figure('Name', 'Figure 3.2: Spatial Phase Evolution', 'Color', 'w', 'Position', [100, 100, 900, 600]);
hold on; grid on;
runs_to_plot = [9, 8, 7, 6, 5, 4, 3, 2, 1]; % Selected runs across the full spectrum
% --- Generate maximally separated colours ---
n = length(runs_to_plot);
base_colors = lines(n);

order = [];
left = 1;
right = n;

while left <= right
    order = [order left];
    if left < right
        order = [order right];
    end
    left = left + 1;
    right = right - 1;
end

colors_matrix = [
    0.0000 0.4470 0.7410  % blue
    0.8500 0.3250 0.0980  % orange
    0.9290 0.6940 0.1250  % yellow
    0.4940 0.1840 0.5560  % purple
    0.4660 0.6740 0.1880  % green
    0.3010 0.7450 0.9330  % cyan
    0.6350 0.0780 0.1840  % dark red
    0.0000 0.0000 0.0000  % black
    0.7500 0.7500 0.7500  % grey
];

for idx = 1:length(runs_to_plot)
    pp = runs_to_plot(idx);
    label_text = sprintf('Slope %d: %.3f m/km', pp, Slopes(pp) * 10);
    
    plot(AllBatchData.(showcase_batch)(pp).x_km, ...
     AllBatchData.(showcase_batch)(pp).Delft3D.kimag * 1e5, ...
     'Color', colors_matrix(idx,:), ...
     'LineWidth', 2, ...
     'DisplayName', label_text);
end
xlim([0, 400]); 
ylim([0.7, 2.5]);
xlabel('Distance along estuary (km)', 'FontWeight', 'bold'); 
ylabel('k_{imag} (\times10^{-5} rad/m) (phase shift)', 'FontWeight', 'bold');
title('Evolution of k_{imag} across Transition Steepness', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
saveas(fig32, fullfile(saveDirMain, 'Figure_3_2_kimag_Evolution.png'));

%% --- FIGURE 3.3: Spatial Damping Evolution (kreal Profiles) ---
fig33 = figure('Name', 'Figure 3.3: Spatial Damping Evolution', 'Color', 'w', 'Position', [120, 120, 900, 600]);
hold on; grid on;

for idx = 1:length(runs_to_plot)
    pp = runs_to_plot(idx);
    label_text = sprintf('Slope %d: %.3f m/km', pp, Slopes(pp) * 10);
    
    plot(AllBatchData.(showcase_batch)(pp).x_km, ...
     AllBatchData.(showcase_batch)(pp).Delft3D.kreal * 1e5, ...
     'Color', colors_matrix(idx,:), ...
     'LineWidth', 2, ...
     'DisplayName', label_text);
end
xlim([0, 400]); 
xlabel('Distance along estuary (km)', 'FontWeight', 'bold'); 
ylabel('k_{real} (\times10^{-5} rad/m) (amplitude damping) ', 'FontWeight', 'bold');
title('Evolution of k_{real} across Transition Steepness', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
saveas(fig33, fullfile(saveDirMain, 'Figure_3_3_kreal_Evolution.png'));

%% --- FIGURE 3.4: Analytical Misfit & Reflection Threshold ---
fig34 = figure('Name', 'Figure 3.4: Reflection Threshold Metric', 'Color', 'w', 'Position', [150, 150, 900, 600]);
hold on; grid on;

marker_styles = {'o-', 's-', 'd-'};
curve_colors  = {[0.85 0.33 0.1], [0.0 0.45 0.74], [0.47 0.67 0.19]}; % Orange, Blue, Green


% 1. Plot the multi-batch curves
for b = 1:length(batches)
    plot(Slopes, mean_error_seaward(b, :), marker_styles{b}, 'Color', curve_colors{b}, ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', curve_colors{b}, ...
         'DisplayName', sprintf('Chézy = %d', chezy_values(b)));
end

% 2. Clean up axes and layout styling
set(gca, 'XDir', 'reverse', 'XScale', 'log', 'Box', 'off'); 
xlabel('Bed Slope (m/km)', 'FontWeight', 'bold');
ylabel('Mean Seaward Misfit Error \epsilon_k (%) [0-150 km]', 'FontWeight', 'bold');
title('Analytical Vs. Modelled Discrepancies with Bed Steepness', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');

% 3. Add normalized directional indicators
text(0.02, 0.05, 'Steep slope', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
text(0.98, 0.05, 'Gentle slope', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'right');

% 4. Dynamic text label positioning for Run numbers (anchored below the Chézy 30 curve)
chezy30_idx = find(chezy_values == 30); 
if isempty(chezy30_idx); chezy30_idx = 1; end % Fallback to first curve if index match fails
max_err = max(mean_error_seaward(:));

for pp = 1:length(Slopes)
    % Subtracted offset and forced Top alignment to sit neatly below the red line
    text(Slopes(pp), mean_error_seaward(chezy30_idx, pp) - max_err * 0.06, sprintf('Slope %d', pp), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
         'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
end

drawnow;
saveas(fig34, fullfile(saveDirMain, 'Figure_3_4_Reflection_Threshold.png'));

fprintf('\nProcessing complete! All 4 structured figures successfully generated and saved inside:\n -> %s\n', saveDirMain);

%% --- FIGURE 3.6: Reconstructed Tidal Amplitude Profiles (0 to 700 km) ---
fig36 = figure('Name', 'Figure 3.6: Intuitive Amplitude Comparison', 'Color', 'w', 'Position', [100, 100, 1100, 850]);

target_runs =[9, 5, 4, 1]; % Runs corresponding to Gentle, Intermediate, Sharp, Steep
run_labels = {'Slope 9 (Gentle Slope)', 'Slope 5 (Intermediate Slope)', 'Slope 4 (Sharp Slope)', 'Slope 1 (Steep Step)'};
slope_texts = {'0.02 m/km', '0.1137 m/km', '0.20 m/km', '1.00 m/km'};

% Boundary condition matching: your seaward boundary forcing amplitude is h0 = 1.0m
A0 = 1.0; 

for col = 1:4
    pp = target_runs(col);
    
    % Extract physical spatial grid data for this specific run
    x_km = AllBatchData.(showcase_batch)(pp).x_km;
    x_m  = x_km * 1000; % Convert km to meters for correct integration step scale
    
    % Crop indices to isolate the domain up to the 700 km mark (before damping reach)
    idx_700 = x_km <= 700;
    x_plot = x_km(idx_700);
    x_m_int = x_m(idx_700);
    
    % Extract analytical k_real curves (amplitude damping factors)
    kr_k0  = AllBatchData.(showcase_batch)(pp).k0.kreal(idx_700);
    kr_kf  = AllBatchData.(showcase_batch)(pp).kf.kreal(idx_700);
    kr_kj  = AllBatchData.(showcase_batch)(pp).kj.kreal(idx_700);
    
    % Reconstruct analytical amplitudes via cumulative spatial integration 
    amp_k0  = A0 * exp(-cumtrapz(x_m_int, kr_k0));
    amp_kf  = A0 * exp(-cumtrapz(x_m_int, kr_kf));
    amp_kj  = A0 * exp(-cumtrapz(x_m_int, kr_kj));
    
    % Extract actual Delft3D amplitude directly from the saved k_real tracking
    kr_d3d = AllBatchData.(showcase_batch)(pp).Delft3D.kreal(idx_700);
    amp_d3d = A0 * exp(-cumtrapz(x_m_int, kr_d3d));
    
    % Subplot matrix layout (2 rows x 2 columns)
    subplot(2, 2, col); hold on; grid on;
    
    % Plot curves with matching color theme from Figure 3.1
    plot(x_plot, amp_d3d, 'b', 'LineWidth', 2);
    plot(x_plot, amp_k0,  'g--');
    plot(x_plot, amp_kf,  'r-.');
    plot(x_plot, amp_kj,  'm:');
    
    % Draw a visual indicator at x = 150 km where the step begins
    xline(150, 'k:', 'Onset of Slope', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'top');
    
    % Formatting Layout
    xlim([0, 700]);
    ylim([0, 1.2]); % Allows visual space to see partial reflection pooling before the step
    xlabel('Distance along estuary x (km)', 'FontWeight', 'bold');
    title(run_labels{col}, 'FontSize', 11, 'FontWeight', 'bold');
    
   % Text Box for Slope Values
    text(0.95, 0.92, ['Slope: ' slope_texts{col}], ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'FontSize', 9, 'FontAngle', 'italic', 'BackgroundColor', 'w');
    
    if col == 1 || col == 3
        ylabel('Tidal Amplitude A (m)', 'FontWeight', 'bold');
    end
    
    if col == 1
        legend('Delft3D', 'k_0 (Frictionless)', 'k_f (Frictional)', 'k_j (Jay)', ...
               'Location', 'east');
    end
end

% Save out to the main-text figures folder
saveas(fig36, fullfile(saveDirMain, 'Figure_3_6_Physical_Amplitude_Profiles.png'));

%%
fig37 = figure('Name', 'Figure 3.7: Spatial Phase Evolution (Radians)', 'Color', 'w', 'Position', [100, 100, 950, 550]);
hold on; grid on;

runs_to_plot =[9,8,7,6,5,4,3,2,1]; 
showcase_batch = 'Final_friction_45'; 

colors_matrix = [
    0.0000 0.4470 0.7410;  % blue
    0.8500 0.3250 0.0980;  % orange
    0.9290 0.6940 0.1250;  % yellow
    0.4940 0.1840 0.5560;  % purple
    0.4660 0.6740 0.1880;  % green
    0.3010 0.7450 0.9330;  % cyan
    0.6350 0.0780 0.1840;  % dark red
    0.0000 0.0000 0.0000;  % black
    0.7500 0.7500 0.7500   % grey
];

for idx = 1:length(runs_to_plot)
    pp = runs_to_plot(idx);
    
    % Extract spatial components
    x_km = AllBatchData.(showcase_batch)(pp).x_km;
    x_m  = x_km * 1000; % Distance in meters for correct integration scaling
    kimag_rad_m = AllBatchData.(showcase_batch)(pp).Delft3D.kimag;
    
    % Numerically integrate k_imag (rad/m) over distance (m) to get total radians
    phase_radians = cumtrapz(x_m, kimag_rad_m);
    
    label_text = sprintf('Slope %d: %.3f m/km', pp, Slopes(pp) * 10);
    
    % Plot the integrated phase profile in radians
    plot(x_km, phase_radians, ...
         'Color', colors_matrix(idx,:), ...
         'LineWidth', 2, ...
         'DisplayName', label_text);
end

% Layout modifications for Radians
xlim([0, 275]); 
ylim([0, 6]); % Fits the cumulative radian accumulation from 0 to 400 km
xlabel('Distance along estuary (km)', 'FontWeight', 'bold'); 
ylabel('Phase \phi (rad)', 'FontWeight', 'bold');
title('S_2 Spatial Phase Evolution (Chézy = 45)', 'FontSize', 12, 'FontWeight', 'bold');

% Draw a visual guideline for the transition start zone
xline(150, 'k:', 'Transition Start', ...
      'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');

% Positions legend cleanly outside on the right edge (or northwest as requested)
legend('Location', 'northwest'); 

saveas(fig37, fullfile(saveDirMain, 'Figure_3_7_Phase_Shift_Radians.png'));

%%
%% --- FIGURE 3.8: Multi-Slope Time Series at Key Stations ---
fig38 = figure('Name', 'Figure 3.8.1: Tidal Evolution Across Slopes', 'Color', 'w', 'Position', [100, 100, 1200, 800]);

target_runs =[9,5,4,1]; % Representative selection spanning gentle to steep
stations_km =[0,50,100,150,200,250]; 
showcase_batch = 'Final_friction_45';

colors_matrix = [
    0.0000 0.4470 0.7410;  % Slope 9: Blue
    0.8500 0.3250 0.0980;  % Slope 8
    0.9290 0.6940 0.1250;  % Slope 7: Yellow
    0.4940 0.1840 0.5560;  % Slope 6
    0.4660 0.6740 0.1880;  % Slope 5: Green
    0.3010 0.7450 0.9330;  % Slope 4: Cyan
    0.6350 0.0780 0.1840;  % Slope 3
    0.0000 0.0000 0.0000;  % Slope 2
    0.7500 0.7500 0.7500   % Slope 1: Grey
];

baseDir = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\' showcase_batch filesep];

for s = 1:length(stations_km)
    subplot(3, 2, s); hold on; grid on;
    target_x_km = stations_km(s);
    
    for idx = 1:length(target_runs)
        pp = target_runs(idx);
        runName = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', abs(-Qall), h0*100, lcAll, pp);
        matFile = [baseDir runName filesep 'md.mat'];
        
        if ~isfile(matFile), continue; end
        
           simData = load(matFile, 'md');
        
        if idx == 1 && s == 1
           
            time_hours = simData.md.time / 60; 
        end
        
      
        x_m_loc = simData.md.x(:, 2);
        [~, idx_station] = min(abs(x_m_loc - (target_x_km * 1000)));
        
        
        wlv_series = squeeze(simData.md.wlv(:, idx_station, 2));
        
       
        label_text = sprintf('Slope %d: %.3f m/km', pp, Slopes(pp) * 10);
        
        plot(time_hours, wlv_series, ...
             'Color', colors_matrix(pp,:), ...
             'LineWidth', 1.5, ...
             'DisplayName', label_text);
    end
    
    title(sprintf('Station x = %d km', target_x_km), 'FontWeight', 'bold');
    xlim([max(time_hours)-36, max(time_hours)]); % Zoom smoothly on the final 36 hours (1.5 days)
    ylim([-1.2, 1.2]);
    
    if mod(s, 2) ~= 0, ylabel('\eta (m)', 'FontWeight', 'bold'); end
    if s > 4, xlabel('Time (hours)', 'FontWeight', 'bold'); end
    if s == 1, legend('Location', 'northwest', 'FontSize', 8); end
end

saveas(fig38, fullfile(saveDirMain, 'Figure_3_8_1_Station_Time_Series.png'));
%%
%% ========================================================================
%  FIGURE 3.9: Seaward Error Analysis Against Max Delta
%  ========================================================================
fig39 = figure('Name', 'Figure 3.9: Misfit Error Against Delta Parameter', 'Color', 'w', 'Position', [100, 100, 750, 500]);
hold on; grid on;


marker_styles = {'o-', 's-', 'd-'};
batch_legends = {'Chézy = 30', 'Chézy = 45', 'Chézy = 60'};
batch_colors = {[0.85 0.33 0.1], [0.0 0.45 0.74], [0.47 0.67 0.19]}; % Orange, Blue, Green

for b = 1:length(batches)
    batch = batches{b};
    max_delta_per_run = zeros(length(Slopes), 1);
    
    for pp = 1:length(Slopes)
        % Extract the exact spatial grid properties stored for this run
        x_m = AllBatchData.(batch)(pp).x_km * 1000;
        h_local = AllBatchData.(batch)(pp).h_local;
        Nx = length(x_m);
        
        % --- Reconstruct the Exact Analytical Slope Metrics ---
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
        
     
        Delta = zeros(Nx, 1);
        for xx = 1:Nx
            if inv_La(xx) ~= 0
                La_val = 1 / inv_La(xx);
                Delta(xx) = sqrt(g * h_local(xx)) / (2 * La_val * w);
            end
        end
        
        % Extract the peak Delta value inside the active transition zone
        max_delta_per_run(pp) = max(Delta(idx_start:idx_end));
    end
    
    % Plot mean seaward error against the true structural peak Delta
    plot(max_delta_per_run, mean_error_seaward(b, :), ...
         marker_styles{b}, 'Color', batch_colors{b}, ...
         'LineWidth', 2, 'MarkerFaceColor', batch_colors{b}, ...
         'MarkerSize', 7, 'DisplayName', batch_legends{b});
end

% Layout Adjustments
xlabel('Topographic Convergence Parameter \Delta_{max}', 'FontWeight', 'bold');
ylabel('Mean Seaward Misfit Error \epsilon_k (%)', 'FontWeight', 'bold');
title('Analytical Misfit Scaling with Peak Convergence Parameter (\Delta)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest');

saveas(fig39, fullfile(saveDirMain, 'Figure_3_9_Error_vs_Delta.png'));


%% ========================================================================
%  4. HELPER CONDITIONAL VALUE ASSIGNMENT FUNCTION
%  ========================================================================
function val = idxchoose(condition, trueVal, falseVal)
    val = zeros(size(condition));
    val(condition) = trueVal(condition);
    if isscalar(falseVal)
        val(~condition) = falseVal;
    else
        val(~condition) = falseVal(~condition);
    end
end
