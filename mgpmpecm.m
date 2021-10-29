%% Function to calculate optimal control, final state, and state vector

function [u_opt, x, lambda] = mgpmpecm(outage_len, pv_outage, ld_outage, ...
    x_max, x_min, inv_rating, inv_eff, ...
    socdot, P_batt_range_socdot, SOC_range_socdot, ...
    dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
    h, lambda_init)
% Function to compute optimal control strategy based on Pontryagin
% Minimization Principle (PMP)
% Inputs:
% outage_len = length of outage period (e.g. one day) in 15 min. intervals
% pv_outage = pv profile during outage period, as kW
% ld_outage = load profile during outage period, as kW
% x_max = upper bound of state
% x_min = lower bound of state
% inv_rating = kWAC rating of BESS inverter
% inv_eff = AC/DC conversion efficiency of BESS inverter
% socdot = map from (pack power, SOC) --> change in SOC
% P_batt_range_socdot = corresponding pack power range for socdot map
% SOC_range_socdot = corresponding SOC range for socdot map
% dfdx = map from (pack power, SOC) --> dSOC_dot/dSOC
% P_batt_range_dfdx = corresponding pack power range for dfdx map
% SOC_range_dfdx = corresponding SOC range for dfdx map
% lambda = initial costate
% h = fraction of the hour
% Outputs: 
% u_opt = optimal control
% x = state as a function of time from t=1:day_len+1
% lambda = costate
x = zeros(outage_len+1, 1);
u_opt = zeros(outage_len, 1);
lambda = zeros(outage_len+1, 1);

% TODO: Change this assumption
x(1) = x_max;

lambda(1) = lambda_init;

% Set penalty (K) for violating constraints
K = 200; % TODO: Sensitivity analysis on this

    for t = 1:outage_len
        % First, work in all AC power
        % Determine control space at time t based on PV, load power
        % assuming that ESS can only charge off of PV, discharge to load
        % But bounded by inverter rating
        u_min_t = -min(pv_outage(t), inv_rating);
        u_max_t = min(ld_outage(t), inv_rating);

        u_range_t = sort([linspace(u_min_t,u_max_t) 0]); % include 0 power as an option

        % Calculate optimal control with control space at time t
        H_t = zeros(length(u_range_t), 1);

        for v = 1:length(u_range_t)
            % Put rules (~ environment) for each dispatch here

            % Initialize PV, load, and ESS energy resources; DG power, penalty
            pv_t = pv_outage(t);
            l_t = ld_outage(t);
            x_t = x(t);
            dg_t = 0;

            u_t = u_range_t(v);

            % First, charge PV off of ESS if u is negative
            pv_t = pv_t - (-min(0,u_t));

            % Dispatch portion of load from ESS if u is positive
            l_t = l_t - (max(0,u_t));

            % Balance the load and PV, with excess PV curtailed or excess load
            % supplied by DG
            if pv_t > l_t
                % Subtract off load from PV and curtail the rest
                % DG is not used here
            elseif pv_t < l_t
                % Subtract off PV from load and supply remainder load with DG
                l_t = l_t - pv_t;
                dg_t = l_t;
            else
                % Load and PV are in perfect balance! How??
                % DG is not used here
            end

            % Calculate penalty contribution
            % First, convert AC power to DC power
            if u_t >= 0
                u_DC_t = u_t/inv_eff;
            else
                u_DC_t = u_t * inv_eff;
            end
            
            % Use socdot map to find SOC at next time step t+1
            socdot_t = interp2(P_batt_range_socdot, SOC_range_socdot, socdot, u_DC_t, x_t, 'spline');
            
            % Find SOC at next time step t+1
            x_t1 = x_t + socdot_t*h;
            
            % Calculate penalty
            if x_t1 > x_max
                w = K;
            elseif x_t1 < x_min
                w = -K;
            elseif x_min <= x_t1 && x_t1 <= x_max
                w = 0;
%             else
%                 w = K;
            end
            
            % Calculate fuel consumption from function
            fuel_cons_L = genset_model(dg_t); 

            % Put resulting Hamiltonian here
            H_t(v) = fuel_cons_L + (lambda(t) + w) * socdot_t;
        end
        
        % Get index of H_t where it is minimized;
        [~, min_ind] = min(H_t);

        % Find the corresponding BESS AC power value as the optimal control
        u_opt_t = u_range_t(min_ind);
        
        % Store optimal BESS AC power control
        u_opt(t) = u_opt_t;

        % First, convert AC power to DC power
        if u_opt_t >= 0
            u_DC_opt_t = u_opt_t/inv_eff;
        else
            u_DC_opt_t = u_opt_t * inv_eff;
        end

        % Use socdot map to find SOC at time step t+1
        socdot_t = interp2(P_batt_range_socdot, SOC_range_socdot, socdot, u_DC_opt_t, x(t), 'spline');
        
        % Update the SOC for the next time step
        x(t+1) = x(t) + socdot_t*h;
        
        % Update the costate
        % 1) Calculate penalty
        if x(t) > x_max
            w = K;
        elseif x(t) < x_min
            w = -K;
        elseif x_min <= x(t) && x(t) <= x_max
            w = 0;
%             else
%                 w = K;
        end
        
        % 2) Compute dfdx
        dfdx_t = -interp2(P_batt_range_dfdx, SOC_range_dfdx, dfdx, u_DC_opt_t, x(t), 'spline');
        
        % 3) Compute lambda_dot
        lambda_dot = -(lambda(t) + w)*dfdx_t;

        % 4) Update the costate for the next time step
        lambda(t+1) = lambda(t) + lambda_dot*h;
    end

    
    % function end
end
