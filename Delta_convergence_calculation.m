% Quick Delta Explorer Tool
g = 9.81; w = 2*pi/(12*3600); h_start = 6.0;

% Type any target bed slope in m/km here to test it
target_slope_m_km = 0.075; 

% Mathematical Conversion
s_m_m = target_slope_m_km / 1000;
Lh_km = (h_start / s_m_m) / 1000;
Delta = (s_m_m * sqrt(g)) / (2 * w * sqrt(h_start));

fprintf('Slope: %.4f m/km | Parameter: %.4f | Lh: %.1f km | Delta: %.4f | 2*Delta: %.4f\n', ...
    target_slope_m_km, target_slope_m_km/10, Lh_km, Delta, 2*Delta);
