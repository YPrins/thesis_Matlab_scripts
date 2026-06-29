function [amp1,amp2,phase1,phase2]= get_S2_S4(signal,t,w)

    M = [cos(w.*t') sin(w.*t') cos(2*w.*t') sin(2*w.*t')];
   
    c = M \ signal;  

   % Extract coefficients
    A1 = c(1);  % cos(wt)
    B1 = c(2);  % sin(wt)
    A2 = c(3);  % cos(2wt)
    B2 = c(4);  % sin(2wt)
    % Amplitudes
    amp1 = sqrt(A1^2 + B1^2);  % First harmonic (omega)
    amp2 = sqrt(A2^2 + B2^2);  % Second harmonic (2*omega)
    %amp2 =[];
    % Phases (in radians)
    phase1 = atan2(B1, A1);  % atan2(sin, cos)
    phase2 = atan2(B2, A2);
    %phase2 = [];


end 