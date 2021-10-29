%% Function to compute derivative of cell current wrt internal resistance
% Kevin Moy
% 9/1/21
% Inputs: 
% x, SOC single value between [0.095, 1]
% u, battery DC power in kilowatts
% s, number of cells in series
% p, number of cells in parallel

function dIcell_dR0 = didr(x, u, s, p)
    % Compute V_oc, R_0 from function approximations
    v = voc(x);
    r = r0(x);
    
    % convert to cell power IN WATTS
    u_cell = u/(s*p)*1000;
    
%     sqrt_term = sqrt( ( (v)/(2*r) )^2 - (u_cell)/(r) );
    sqrt_term = sqrt(v^2 - 4*r*u_cell);
    
%     dIcell_dR0 = -(v)/(2*r^2) - ...
%         ( ( (u_cell)/(r^2) ) - ( (v^2)/(2*r^3) ) ) / (2 * sqrt_term);

    dIcell_dR0 = u_cell/(r*sqrt_term) - (v-sqrt_term)/(2*r^2);
end