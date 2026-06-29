function [depth, cor, wlv, MN] = make_depth_trapjes_stairs(ls, dx, dy, bc, hd, lc, slope, num_steps)
%% OFFshore grid
Mdir = bc/dy;
Ndir = (lc(2)/dx +(ls-lc(2))/(2*dx)) +1;

%% Make mesh  
x       = [0:dx:lc(2) lc(2)+dx*2:dx*2:ls];
Xarray  = x.*ones(Mdir+1,1);

Yarray  = nan(Mdir+1,Ndir);
Yarray(1:end,:) = repmat([bc:-dy:0]',1,size(Yarray,2));
Xarray(isnan(Yarray)) = nan; 

%% Save coordinates
cor.x =  Xarray;
cor.y =  Yarray;

%% Channel depth
dD    = (hd(2)-hd(1));          
dL    = dD./(slope/100);        
hd1_id = floor(lc(1)/dx);
hd2_id = ceil((lc(1)-dL)/dx);

depth                    = nan(size(Xarray));
depth(:,1)               = hd(1);
depth(:,hd1_id)          = hd(1);
depth(:,hd2_id)          = hd(2);
depth(:,floor(lc(2)/dx)) = hd(2);
depth(:,end)             = hd(3);
depth = fillmissing(depth,"linear",2);

%% ========================================================================
%  HORIZONTAL STEP DISCRETIZATION
%  =======================================================================
if num_steps > 0
    z_smooth = depth(1, :);
    
    % Identify exact bounding indices of the bathymetric transition slope
    idx_start = min(hd1_id, hd2_id);
    idx_end   = max(hd1_id, hd2_id);
    num_elements = idx_end - idx_start + 1;
    
    if num_elements > num_steps
        % Divide the horizontal coordinate stretch into equal distance chunks
        step_chunk = num_elements / num_steps;
        
        for s = 1:num_steps
            % Map out the exact index span for this step tread block
            i_start = round(idx_start + (s-1) * step_chunk);
            i_end   = round(idx_start + s * step_chunk) - 1;
            if s == num_steps, i_end = idx_end; end % Guard end boundary
            
            % Locate the exact horizontal center of this specific step block
            i_mid = round((i_start + i_end) / 2);
            
            % Sample the smooth depth at the center point
            % This guarantees that the step volume perfectly averages out to the slope
            midpoint_depth = z_smooth(i_mid);
            
            % Flatten the entire distance block into a flat tread layer
            z_smooth(i_start:i_end) = midpoint_depth;
        end
    end
    
    % Re-expand the modified profile back to populate the full 2D mesh grid
    depth = repmat(z_smooth, size(depth, 1), 1);
end
% ========================================================================

depth(isnan(Yarray)) = -999; 

%% ini waterlevel
wlv = repmat(max(1.5, -depth(1,:)+.25), Mdir+1, 1);
wlv(isnan(Xarray)) = -999; 

%% MN
MN = [size(Xarray)];

end