clc;
clear all;

%% ============================================================
%  FOLDERS
%  ============================================================
workfolder  = 'C:\Users\ymarp\OneDrive\Documenten\Thesis\MatLab to use';

batch       = 'Test_stairs_8';

path        = ['C:\Users\ymarp\OneDrive\Documenten\Thesis\Runs\' batch filesep];

%Change line 146 to number of steps wanted
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
Qall = 100;      % discharge values
Hall = [1.00];            % tidal amplitude
S    = [30,0];            % salinity boundary values


%% ============================================================
%  DOMAIN SIZE
%  ============================================================
ls = 1200000;   % river length (m)(1200 km)
bc = 500;      % river width (m)

dx = 100;      % grid resolution x-direction
dy = 500;      % grid resolution y-direction


%% ============================================================
%  BATHYMETRY
%  ============================================================
hd = [12,6,-2];     % depths (Hb1, Hb2, Hb3)
lcAll  = [150000];             % characteristic length
Slopes = [0.1,0.05,0.02,0.01137,0.005,0.002];  % bed slopes
% Keep slopes 1, 2, 4, 5, 8, 9 of original slopes

%% ============================================================
%  PHYSICAL PARAMETERS
%  ============================================================
siglay = 1;      % number of vertical layers
Cin    = 45;      % Chezy friction coefficient

figure; hold on;
%% ============================================================
%  PARAMETER LOOPS
%  ============================================================
for pp = 1:length(Slopes)

for hh = 1:length(lcAll)

for qq = 1:length(Qall)

    %% --------------------------------------------------------
    %  SET RUN PARAMETERS
    %  --------------------------------------------------------
    
    lc = [lcAll(hh),700000];

    cd(workfolder)

    Q1 = -Qall(qq);       % river discharge
    h0 = Hall(1);         % tidal amplitude


    %% --------------------------------------------------------
    %  RUN NAME
    %  --------------------------------------------------------
    
    clear Bct Bcc Bnd

    run = sprintf('run_Qr%d_h0%d_zb1%d_slope%d', ...
        abs(Q1), h0*100, lcAll(hh), pp);


    %% --------------------------------------------------------
    %  FILE PATHS
    %  --------------------------------------------------------
    
    npath = [path filesep run filesep];
    mkdir(npath);

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
    %  BOUNDARY CONDITIONS
    %  ========================================================
    
    bndname = {'SeaBound','RiverBound'};

    % Boundary type
    % Z = water level
    % T = total discharge
    bndtype = {'Z','T'};

    datatype = {'T','T'};    % timeseries

    bndalfa = [0,0,0,0];

    bndvertp = {'Uniform','Uniform'};
    bndvertp_bcc = 'Uniform';

    bndlabA = {[],[]};
    bndlabB = {[],[]};

    % tidal constituents
    cnstit1 = {'M2'};
    A1      = [1];
    phi1    = [0];

    C      = Cin;     % Chezy roughness
    Vicouv = 15;

    t0 = 5760;       % start of mapfile output


    %% ========================================================
    %  CREATE GRID AND DEPTH
    %  ========================================================
    
    use_steps = true;       % Set to true to run stair-step configuration, false for smooth
    num_steps = 8;          % Adjust the number of discrete steps 
    
    if use_steps
        [depth, cor, wlv, MN] = make_depth_trapjes_stairs( ...
            ls, dx, dy, bc, hd, lc, Slopes(pp), num_steps);
        
       
        [depth_smooth, ~, ~, ~] = make_depth_trapjes( ...
            ls, dx, dy, bc, hd, lc, Slopes(pp));
    else
       
        [depth, cor, wlv, MN] = make_depth_trapjes( ...
            ls, dx, dy, bc, hd, lc, Slopes(pp));
    end

    fig = figure('Visible','off', 'Color', 'w');

    % Plot the active bathymetry matrix column layout
    plot(cor.x(1,:)/1000, depth(1,:), ...
        'Color',[0.76 0.60 0.42], ...
        'LineWidth',2, ...
        'DisplayName', 'Active Bathymetry')
    
    hold on
    grid on

    % If running stairs, overlay the smooth baseline line underneath to verify the midpoint matching
    if use_steps
        plot(cor.x(1,:)/1000, depth_smooth(1,:), ...
            'k--', 'LineWidth', 1.2, 'DisplayName', 'Continuous Slope Trend')
    end

    xlabel('Distance along estuary (km)')
    ylabel('Bed level (m)')

    % Sea and river boundary labels
    text(0, max(depth(:)), '  Sea boundary', ...
        'Color','b', ...
        'VerticalAlignment','bottom');

    text(ls/1000 - 45, max(depth(:)), 'River boundary  ', ...
        'Color','b', ...
        'VerticalAlignment','bottom');

     if use_steps
        title(sprintf('Bathymetry Matrix Comparison (Slope = %.4f, Steps = %d)', ...
            Slopes(pp), num_steps), 'FontWeight', 'bold')
        img_name = sprintf('bathymetry_steps_%d', num_steps);
    else
        title(sprintf('Bathymetry Profile (Smooth Slope = %.4f, l_c = %d m)', ...
            Slopes(pp), lcAll(hh)), 'FontWeight', 'bold')
        img_name = 'bathymetry_smooth';
    end

    legend('Location','best', 'Box', 'off')
    set(gca,'YDir','reverse', 'Box', 'off')

       saveas(fig, fullfile(npath, [img_name '.png']))
    print(fig, fullfile(npath, img_name), '-dpng', '-r300')

    close(fig)
    %% ========================================================
    %  WRITE BOUNDARY STRUCTURE
    %  ========================================================
    
    dummy  = 1;
    Mcells = MN(1) + dummy;
    Ncells = MN(2) + dummy;

    McellsC = (bc)/dy + dummy;

    bndnm = {
        [Mcells-1;1;Mcells-McellsC+1;1], ...
        [Mcells-1;Ncells;Mcells-McellsC+1;Ncells]
        };

    Bnd = makebndstruct( ...
        BcName, bndnm, bndname, bndtype, datatype, ...
        bndalfa, bndvertp, bndlabA, bndlabB);


    %% ========================================================
    %  TIDAL FORCING
    %  ========================================================
    
    TimeIn = 0:dt:EndTime;

    w = 2*pi/T;

    EtaIn = -h0 .* cos(w .* TimeIn);   % water level
    QIn   = Q1;                        % discharge


    %% ========================================================
    %  CREATE BOUNDARY FILES
    %  ========================================================
    
    Bct = makebctstruct_tidalseries( ...
        BctName, bndvertp, bndname, ...
        RefTime, TimeUnit, TimeIn, ...
        EtaIn, QIn, -QIn);

    Bcc = makebccstruct( ...
        BccName, bndvertp_bcc, bndname, ...
        RefTime, TimeUnit, EndTime, [S]);


    %% ========================================================
    %  INITIAL SALINITY FIELD
    %  ========================================================
    
    Lmax       = 80000;
    NcellsSmax = floor(Lmax ./ dx);

    Sx = linspace(30,0,NcellsSmax);

    Data.wlv = ones([Mcells,Ncells]);
    Data.wlv(2:end,1:end-1) = wlv;

    Data.u = zeros([Mcells,Ncells,siglay]);
    Data.v = zeros([Mcells,Ncells,siglay]);


    %% ========================================================
    %  WRITE MODEL INPUT FILES
    %  ========================================================
    
    delft3d_io_grd('write',[npath filesep GridName],cor.x,cor.y);

    movefile([npath,GridName], ...
             [npath,GridName,'.enc'])

    delft3d_io_dep('write', ...
        [npath filesep DepName], depth', 'location','cor');

    delft3d_io_bnd('write',[npath filesep BcName],Bnd);

    bct_io('write',[npath,BctName],Bct);
    %bct_io('write',[npath,BccName],Bcc);

    delft3d_io_ini('write',[npath filesep IniName],Data)

    G = delft3d_io_grd('read',[npath filesep GridName '.grd']);


    %% ========================================================
    %  RUNID FILE
    %  ========================================================
    
    fileID = fopen([npath filesep 'runid.'],'w');
    fprintf(fileID,run);
    fclose(fileID);


    %% ========================================================
    %  MASTER DEFINITION FILE (.MDF)
    %  ========================================================
    
    Mdf = delft3d_io_mdf('new','template_gui.mdf');

    Mdf.keywords.ident  = 'Delft3D-FLOW 3.59.01.48550';

    Mdf.keywords.filcco = [GridName,'.grd'];
    Mdf.keywords.filgrd = [GridName,'.enc'];

    Mdf.keywords.mnkmax = [G.mmax,G.nmax,siglay];

    Mdf.keywords.thick  = repmat(1/siglay*100,siglay,1);

    Mdf.keywords.fildep = DepName;

    Mdf.keywords.itdate = '2000-01-01';

    Mdf.keywords.tunit  = 'M';
    Mdf.keywords.tstart = 0;
    Mdf.keywords.tstop  = EndTime;
    Mdf.keywords.dt     = 1;

    Mdf.keywords.zeta0  = 0;

    %Mdf.keywords.sub1   = '   ';
    %Mdf.keywords.sub2   = '   ';

    Mdf.keywords.filbnd = BcName;
    Mdf.keywords.filbct = BctName;
    %Mdf.keywords.filbcc = BccName;

    Mdf.keywords.filic  = IniName;

    Mdf.keywords.roumet = 'C';

    Mdf.keywords.ccofu  = C;
    Mdf.keywords.ccofv  = C;

    Mdf.keywords.rettis = '';
    Mdf.keywords.rettib = '';

    Mdf.keywords.vicouv = 1;
    Mdf.keywords.dicouv = 10;

    Mdf.keywords.vicoww = 0.0001;
    Mdf.keywords.dicoww = 0.0001;

    Mdf.keywords.tkemod = 'K-epsilon';

    Mdf.keywords.phhydr = 'YYYYYY';
    Mdf.keywords.phderv = 'YYY';
    Mdf.keywords.phproc = 'YYYYYYYYYY';
    Mdf.keywords.phflux = 'YYYY';

    Mdf.keywords.flmap = [EndTime*(35/40),5,EndTime];
    Mdf.keywords.flhis = [EndTime*(35/40),5,EndTime];

    Mdf.keywords.flrst = 0;

    delft3d_io_mdf('write',[npath,MdfName],Mdf.keywords);

    fclose('all');

    plot(cor.x,depth(1,:))

end
end
end

