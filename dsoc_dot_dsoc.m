function [dfdx, P_batt_range, SOC_range] = dsoc_dot_dsoc(max_P, min_P, max_SOC, min_SOC, s, p, Q_nom)
% Function to generate map (dSOC_dot/dSOC)

%TODO: make these an argument of the function, but for now leave these
%points as default
numpts_soc = 101;
numpts_P = 51;

% Range of P_batt:
P_batt_range = linspace(min_P,max_P,numpts_P);

% Range of SOC:
SOC_range = linspace(min_SOC,max_SOC,numpts_soc);

% [X, Y] =  meshgrid(P_batt_range, SOC_range)

    for i = 1:numpts_soc
        soc = SOC_range(i);
        for j = 1:numpts_P
            P_batt = P_batt_range(j);
            di_dvoc = didv(soc, P_batt, s, p);
            di_dr0 = didr(soc, P_batt, s, p);
            [~, dvoc_dsoc] = voc(soc);
            [~, dr0_dsoc] = r0(soc);

            dfdx(i,j) = -(1/Q_nom) * (di_dvoc*dvoc_dsoc + di_dr0*dr0_dsoc);
        end
    end

end