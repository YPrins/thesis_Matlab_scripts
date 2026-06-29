function [depth, cor,wlv,MN] = make_depth_trapjes(ls,dx,dy,bc,hd,lc,slope)

%% OFFshore grid (exponential increasing gridsize)
Mchannel = bc/dy;
Mdir     = bc/dy;
Ndir     = (lc(2)/dx +(ls-lc(2))/(2*dx)) +1;
%% Make mesh  
x       = [0:dx:lc(2) lc(2)+dx*2:dx*2:ls];
Xarray  = x.*ones(Mdir+1,1);



Yarray  = nan(Mdir+1,Ndir);
Yarray(1:end,:) = repmat([bc:-dy:0]',1,size(Yarray,2));
Xarray(isnan(Yarray)) = nan; 
% 
% 
% dyyy          = 5000*exp(-x./10000);
% Yarray(1,:) = Yarray(1,:)+ dyyy ;
%% Save cordiantes
cor.x =  Xarray;
cor.y =  Yarray;
%% Channel depth
dD    = (hd(2)-hd(1));          % change in depht
dL    = dD./(slope/100);        % lenght of slope
hd1_id = floor(lc(1)/dx);
hd2_id = ceil((lc(1)-dL)/dx);

depth                    = nan(size(Xarray));
depth(:,1)               = hd(1);
depth(:,hd1_id)          = hd(1);
depth(:,hd2_id)          = hd(2);
depth(:,floor(lc(2)/dx)) = hd(2);
depth(:,end)             = hd(3);
depth = fillmissing(depth,"linear",2);
depth(isnan(Yarray)) = -999; 

%% ini waterlevel
wlv = ones(size(depth));

% takes 0 wlv or 1 meter water level. 
wlv   = repmat(max(1.5, -depth(1,:)+.25),Mdir+1,1);
wlv (isnan(Xarray)) = -999; 
%% MN

MN = [size(Xarray)];
end