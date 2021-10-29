%% Function to compute open-circuit voltage and its derivative wrt SOC
% Kevin Moy
% 9/1/21
% Inputs: 
% SOC, single value between [0.095, 1]
%%
function [v_oc, dv_oc] = voc(soc)
    % TODO: Better error handling
    if soc >= 0.0949 && soc <= 0.1885
        v_oc = 2.6071*soc + 2.9333;
        dv_oc = 2.6071;
    elseif soc > 0.1885 && soc <= 1
        v_oc = 0.9374*soc + 3.2436;
        dv_oc = 0.9374;
    else
        v_oc = Inf;
        dv_oc = Inf;
    end
end