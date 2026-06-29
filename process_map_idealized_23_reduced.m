clc; clear; close all;

%% ------------------------------------------------------------------------
%  USER SETTINGS (students only edit this block)
% -------------------------------------------------------------------------
symType   = 2;                 % 2 for 2DV models
baseDir   = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\Test_stairs_8\';   % <-- CHANGE THIS
addpath('C:\Users\ymarp\OneDrive\Documenten\Thesis\MatLab to use\process_map_idealized_23_reduced.m\');               % folder containing helper functions
% -------------------------------------------------------------------------

%% Check if base directory exists
if ~isfolder(baseDir)
    error('Base directory not found: %s', baseDir);
end

%% Get list of run folders
fls = dir(baseDir);
fls = fls([fls.isdir]);        % keep only directories

% Remove '.' and '..'
fls = fls(~ismember({fls.name},{'.','..'}));

fprintf('Found %d model runs.\n', numel(fls));

%% Loop through all runs
for ii = 1:numel(fls)

    runName = fls(ii).name;
    fprintf('\nProcessing run: %s\n', runName);

    runDir   = fullfile(baseDir, runName);
    fileTrim = fullfile(runDir, ['trim-' runName '.dat']);
    fileGrd  = fullfile(runDir, [runName '.grd']);

    %% Check files
    if ~isfile(fileTrim)
        warning('Missing trim file for run %s. Skipping...', runName);
        continue
    end
    if ~isfile(fileGrd)
        warning('Missing grid file for run %s. Skipping...', runName);
        continue
    end

    %% Read grid
    Grid = delft3d_io_grd('read', fileGrd);

    %% Read model data
    vsinfo = vs_use(fileTrim, 'quiet');

    %% Time settings
    Tp    = 12 * 60;   % tidal period [min]
    dt    = 5;         % time step [min]
    time  = 0:dt:Tp*5;
    ntuse = numel(time);

    nstps = vsinfo.GrpDat(1).SizeDim;

    if ntuse > nstps
        warning('Requested time window exceeds available model steps. Skipping run %s.', runName);
        continue
    end

    %% Extract model data
    md = get_model_data(vsinfo, nstps-ntuse+1 : nstps);
    md.time = time;

    %% Save output
    save(fullfile(runDir, 'md.mat'), 'md', '-v7.3');

    fprintf('Saved md.mat for run %s\n', runName);

end

fprintf('\nAll runs processed.\n');
