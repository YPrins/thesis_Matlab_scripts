clc;
clear all;
close all;

%% ============================================================
%  USER FOLDER & SYSTEM CONFIGURATION
%  ============================================================
workfolder  = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\MatLab to use';
batch       = 'Stairs_2_visualisation';
base_path   = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs';


path = fullfile(base_path, batch);

% Ensure the base execution output directory exists
if ~exist(path, 'dir')
    mkdir(path);
end


if exist(workfolder, 'dir')
    cd(workfolder);
else
    error('The specified workfolder does not exist: %s', workfolder);
end

%% ============================================================
%  TIME SETTINGS
%  ============================================================
RefTime  = 20000101;      % reference time (minutes since this date)
EndTime  = 60*24*60;      % simulation time (~1 month)
TimeUnit = 'minutes';

dt = 1;                   % timestep
T  = 12*60;               % tidal period


%% ============================================================
%  MODEL PARAMETERS
%  ============================================================
Qall = 100;               % discharge values
Hall = [1.00];            % tidal amplitude


%% ============================================================
%  DOMAIN SIZE
%  ============================================================
ls = 1200000;             % river length (m) (1200 km)
bc = 500;                 % river width (m)

dx = 100;                 % grid resolution x-direction
dy = 500;                 % grid resolution y-direction


%% ============================================================
%  BATHYMETRY SETUP
%  ============================================================
hd     = [12,6,-2];       % depths (Hb1, Hb2, Hb3)
lcAll  =[150000];        % characteristic length
Slopes = [0.002]; % bed slopes


%% ============================================================
%  PHYSICAL PARAMETERS
%  ============================================================
siglay = 1;               % number of vertical layers
Cin    = 45;              % Chezy friction coefficient


%% ============================================================
%  PARAMETER LOOPS
%  ============================================================
for pp = 1:length(Slopes)
    for hh = 1:length(lcAll)
        for qq = 1:length(Qall)

            %% --------------------------------------------------------
            %  SET RUN PARAMETERS
            %  --------------------------------------------------------
            lc = [lcAll(hh), 700000];

            Q1 = -Qall(qq);       % river discharge
            h0 = Hall(1);         % tidal amplitude

            %% --------------------------------------------------------
            %  RUN NAME GENERATION
            %  --------------------------------------------------------
            clear Bct Bcc Bnd
            run = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', abs(Q1), h0*100, lcAll(hh), pp);

            %% --------------------------------------------------------
            %  FILE PATHS AND RUN DIRECTORY CREATION
            %  --------------------------------------------------------
            npath = fullfile(path, run);
            if ~exist(npath, 'dir')
                mkdir(npath);
            end

            BcName   = [run '.bnd'];
            BctName  = [run '.bct'];
            BccName  = [run '.bcc'];
            BcaName  = [run '.bca'];
            DepName  = [run '.dep'];
            IniName  = [run '.ini'];
            ObsName  = [run '.obs'];
            XSName   = [run '.crs'];
            GridName = run;
            MdfName  = [run '.mdf'];

            %% ========================================================
            %  CREATE GRID AND DEPTH
            %  =======================================================
            use_steps = true;       % true = stair-step, false = smooth slope
            num_steps = 2;          % Total target horizontal step intervals

            if use_steps
                [depth, cor, wlv, MN] = make_depth_trapjes_stairs( ...
                    ls, dx, dy, bc, hd, lc, Slopes(pp), num_steps);
                
                [depth_smooth, ~, ~, ~] = make_depth_trapjes( ...
                    ls, dx, dy, bc, hd, lc, Slopes(pp));
            else
                [depth, cor, wlv, MN] = make_depth_trapjes( ...
                    ls, dx, dy, bc, hd, lc, Slopes(pp));
            end

            FONT_NAME = 'Helvetica'; 
            FS_LABEL  = 10;
            FS_TITLE  = 11;
            FS_TEXT   = 9;

            % Setup tight canvas dimension matching a standard single-column (14cm x 10cm)
            fig = figure('Visible', 'off', 'Color', 'w');
            set(fig, 'Units', 'centimeters', 'Position', [2, 2, 14, 10]);

            % Plot 1: Active Bathymetry 
            plot(cor.x(1,:)/1000, depth(1,:), ...
                'Color', [0.85, 0.51, 0.17], ...
                'LineWidth', 1.8, ...
                'DisplayName', 'Stepped Bathymetry');
            hold on;

            % Plot 2: Reference Smooth Curve
            if use_steps
                plot(cor.x(1,:)/1000, depth_smooth(1,:), ...
                    ':', 'Color', [0.29, 0.33, 0.41], ...
                    'LineWidth', 1.4, ...
                    'DisplayName', 'Continuous Bathymetry');
            end

            % Format Canvas Axes & Grid
            grid on;
            set(gca, ...
                'FontName', FONT_NAME, ...
                'FontSize', FS_LABEL-1, ...
                'YDir', 'reverse', ...
                'Box', 'off', ...
                'XMinorTick', 'on', ...
                'YMinorTick', 'on', ...
                'TickDir', 'out', ...
                'LineWidth', 0.8, ...
                'GridLineStyle', ':', ...
                'GridColor', [0.6, 0.6, 0.6], ...
                'GridAlpha', 0.3);

            % Labels & Titles
            xlabel('Distance along estuary (km)', 'FontName', FONT_NAME, 'FontSize', FS_LABEL, 'FontWeight', 'bold');
            ylabel('Bed level (m)', 'FontName', FONT_NAME, 'FontSize', FS_LABEL, 'FontWeight', 'bold');

           
            xLim = xlim;
            yLim = ylim;
            yPos_Text = (yLim(1) + yLim(2)) / 2; 

            text(xLim(1) + (xLim(2)-xLim(1))*0.02, yPos_Text, 'Sea Boundary', ...
                'FontName', FONT_NAME, 'FontSize', FS_TEXT, 'Color', [0.1, 0.3, 0.6], ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

            text(xLim(2) - (xLim(2)-xLim(1))*0.02, yPos_Text, 'River Boundary', ...
                'FontName', FONT_NAME, 'FontSize', FS_TEXT, 'Color', [0.1, 0.3, 0.6], ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

            % Dynamic Titles & Filenames
            if use_steps
                title_str = ('Bathymetry Comparison (Slope = 0.002 (0.02 m/km), Steps = 2)');
                img_name  = sprintf('bathymetry_steps_%d', num_steps);
            else
                title_str = sprintf('Bathymetry Profile (Smooth Slope = %.4f)', Slopes(pp));
                img_name  = 'bathymetry_smooth';
            end
            title(title_str, 'FontName', FONT_NAME, 'FontSize', FS_TITLE, 'FontWeight', 'bold');

         
            lgd = legend('Location', 'best', 'Box', 'off');
            set(lgd, 'FontName', FONT_NAME, 'FontSize', FS_LABEL-1);

          
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left, bottom, ax_width, ax_height];

            % Save Output
            saveas(fig, fullfile(npath, [img_name '.png']));
            print(fig, fullfile(npath, img_name), '-dpng', '-r300');
            print(fig, fullfile(npath, img_name), '-depsc', '-r600'); 

            close(fig);
        end
    end
end
