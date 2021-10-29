%% Function to compute degradation from http://dx.doi.org/10.1016/j.jpowsour.2014.12.047
% (Cordoba et al.)

function [cap_loss] = emp_deg_model(SOC, Q_nom, h, temp_deg_C)
% Inputs:
% SOC = state evolution, unitless
% Q_nom = nominal capacity, in Ah
% h = data resolution, in fraction of hour (unitless)
% temp_deg_C = temperature, in degrees Celsius (average for now)

    % Assume average temperature for now??
    T = temp_deg_C + 273.15; %25 degrees C, in Kelvin

    % Constants:
    alpha_c = 137;
    beta_c = 420;
    gamma_c = 9610;
    b = 0.34;
    c = 3;
    z = 0.48;    
    SOC_0 = 0.25;
    E_a_c = 22406; %J/mol
    R_g = 8.314; %J/(K*mol);
    
    %Threshold to consider charge-sustaining (in A)
    UE_CS = 0.1;
    
    % Compute cell current
    I_cell = -diff(SOC)/h;
    
    % Compute ratio from u
    CD = length(find(I_cell>0));
    CS = length(find(abs(I_cell) < UE_CS));
    Ratio = CD/(CD+CS);
    
    % Compute minimum SOC
    SOC_min = min(SOC);
   
    % Compute throughput
    Ah = sum(abs(I_cell))*Q_nom*h;
    
    % Compute capacity severity factor
    a_c = alpha_c + beta_c*Ratio^b + gamma_c*(SOC_min-SOC_0)^c;
    
    % Compute capacity fade (in Ah)
    cap_loss = a_c*exp(-E_a_c/(R_g*T))*Ah^z;
    
end