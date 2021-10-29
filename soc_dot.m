function [socdot, P_batt_range, SOC_range] = soc_dot(max_P, min_P, max_SOC, min_SOC, s, p, Q_nom)

%TODO: make these an argument of the function, but for now leave these
%points as default
numpts_soc = 101;
numpts_P = 50;

% Range of P_batt:
P_batt_range = [linspace(min_P,max_P,numpts_P) 0]; % Include 0 power
P_batt_range = sort(P_batt_range);

% Range of SOC:
SOC_range = linspace(min_SOC,max_SOC,numpts_soc);

    for i = 1:numpts_soc
        for j = 1:length(P_batt_range)
            P_batt = P_batt_range(j);
            P_cell = P_batt/(s*p)*1000; %CONVERT TO WATTS
            soc = SOC_range(i);
            [v_oc,~] = voc(soc);
            [r_0, ~] = r0(soc);
            socdot(i,j) = -1/Q_nom * ( (v_oc - sqrt(v_oc^2-4*r_0*P_cell))/(2*r_0) ); 
        end
    end
end
