function fuel_cons = genset_model(dg_power)
% Function to produce output of 457kW diesel genset
% Kevin Moy
% Sep 8 2021

% NOTE: Only valid for dg_power in [0, 457]
% Input: Desired DG output electric power, in kWAC 
% Output: Fuel (diesel) consumption, in liters per hour (L/hr)
    fuel_cons = 1.2349e-04 .* dg_power.^2 + 0.1982 .* dg_power + 16.3602;
end