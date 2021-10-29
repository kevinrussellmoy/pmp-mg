%% Function to compute internal resistance and its derivative wrt SOC
% Kevin Moy
% 9/1/21
% Inputs: 
% SOC, single value between [0.095, 1]
%%
function [r_0, dr_0] = r0(soc)
    % TODO: Better error handling
    if soc >= 0.0949 && soc <= 1
        r_0 = 1.2317*soc^6 - 4.3378*soc^5 + 6.1645*soc^4 - 4.5158*soc^3 + ...
            1.805*soc^2 - 0.3779*soc + 0.06567;
        dr_0 = (1.2317*6)*soc^5 - (4.3378*5)*soc^4 + (6.1645*4)*soc^3 - ...
            (4.5158*3)*soc^2 + (1.805*2)*soc - 0.3779;
    else
        r_0 = Inf;
        dr_0 = Inf;
    end
    
%     % constant r0:
%     r_0 = 0.04;
%     dr_0 = 0;
end