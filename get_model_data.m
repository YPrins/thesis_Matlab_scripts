
function md = get_model_data(vsinfo,tsp)
    % x,y, in m,n
    md.nlyrs   = vs_let(vsinfo, 'map-const', 'KMAX',  'quiet'); 
    md.x       = squeeze(vs_let(vsinfo, 'map-const',{0}, 'XZ', {0,0}, 'quiet')); 
    md.y       = squeeze(vs_let(vsinfo, 'map-const',{0}, 'YZ', {0,0}, 'quiet')); 
    
    % u,v,w,
    md.u       = vs_let(vsinfo, 'map-series',{tsp}, 'U1', {0,0,0}, 'quiet');     % horizontal velocity (m)
    md.v       = vs_let(vsinfo, 'map-series',{tsp}, 'V1', {0,0,0}, 'quiet');     % horizontal velocity (n)
    %md.w       = vs_let(vsinfo, 'map-series',{tsp}, 'WPHY', {0,0,0}, 'quiet'); 
    
    % scalar fields 
    %md.rho   = vs_let(vsinfo, 'map-series',{tsp}, 'RHO', {0,0,0}, 'quiet'); 
    %md.s     = vs_let(vsinfo, 'map-series',{tsp}, 'R1', {0,0,0,1}, 'quiet'); 
    %md.Av    = vs_let(vsinfo, 'map-series',{tsp}, 'DICWW', {0,0,0}, 'quiet'); 
    
    % depth, wlv
    md.ah          = squeeze(vs_let(vsinfo, 'map-const',{0}, 'GSQS', {0,0}, 'quiet'));
    md.thck       = vs_let(vsinfo, 'map-const',{0},'THICK', {0}, 'quiet');
    md.wlv        = vs_let(vsinfo, 'map-series',{tsp},'S1', {0,0}, 'quiet');
    md.zb        = squeeze(vs_let(vsinfo, 'map-const',{0},'DPS0', {0,0}, 'quiet'));
end 
