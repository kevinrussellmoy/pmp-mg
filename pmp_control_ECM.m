%% PMP Sizing for Microgrid Energy Storage using Equivalent Circuit Model
% Kevin Moy, 8/15(?)/21
clearvars
close all
clc
set(0,'defaultTextInterpreter','latex');
%% Global vartiables (sizes)
PV_ARRAY_SIZE_KW = 210;   % kWAC rating of the PV array
DIESEL_GEN_SIZE_KW = 457;   % kWAC rating of the diesel generator
LIB_INV_SIZE_KW = 250; % kWAC rating of the LIB (Dynapower) inverter
LIB_INV_MAX_DC_V = 830; % Maximum DC voltage of the LIB inverter
LIB_INV_MIN_DC_V = 530; % Minimum DC voltage of the LIB inverter
LIB_NOM_ENG = 475; % kWh nominal energy of the pack
LIB_EFF_DISCHG = 0.975; % Efficiency of discharging ESS, from Dynapower specs
LIB_EFF_CHG = 0.975; % Efficiency of charging ESS, from Dynapower specs
h = 15/60; % Fraction of the hour

%% LIB Configuration
% State is SOC
% Control is AC power, bounded by kWAC rating of the LIB inverter

% Cell nominal voltage:
v_cell_nom = 3.63; % Volts

% Cell capacity:
Q_nom = 4.85; % Ampere-hours

% Number of cells in series
s = floor(LIB_INV_MAX_DC_V/v_cell_nom); %=207

% Number of cells in parallel
p = floor((LIB_NOM_ENG*1000)/(s*v_cell_nom*Q_nom)); %=130

% State limits
x_min = 0.15;
x_max = 0.85;

% Limits of charging and discharging **DC** power
u_min = -LIB_INV_SIZE_KW*LIB_EFF_CHG;
u_max = LIB_INV_SIZE_KW/LIB_EFF_DISCHG;

% Precompute map of dSOC_dot/dSOC
[dfdx, P_batt_range_dfdx, SOC_range_dfdx] = dsoc_dot_dsoc(u_max, u_min, x_max, x_min, s, p, Q_nom);

% Precompute map of SOC_dot
[socdot, P_batt_range_socdot, SOC_range_socdot] = soc_dot(u_max, u_min, x_max, x_min, s, p, Q_nom);

%% Diesel fuel consumption (use genset_model.m)


%% Datetime vector for plotting
YEAR_START = datetime(2019,1,1,0,0,0);
YEAR_END = datetime(2019,12,31,23,45,0);
yr_dt = (YEAR_START:minutes(15):YEAR_END)';

%% Obtain load and PV data:
load('pv_gen.mat');
load('load_cons.mat');

% Plot load, PV data with datetime x-axis TODO MOVE TO END??
% hFig = figure(1);
% set(hFig, 'Position', [100 100 600 500])
% hold on
% plot(yr_dt, ld)
% plot(yr_dt, pv)

% % Plot residual load per day
% hFig = figure(2);
% set(hFig, 'Position', [100 100 600 500])
% hold on
% plot(yr_dt, ld-pv)

%% Define initial, final time

t_ini_str = '09-Oct-2019 00:00:00';

t_fin_str = '12-Oct-2019 17:45:00';

t_vec = [YEAR_START; datetime([t_ini_str; t_fin_str])];

t_ind = ((minutes(t_vec(2:3) - t_vec(1)))/15) + 1;

%% Select load and PV data

ld1 = ld(t_ind(1):t_ind(2))/2;
pv1 = pv(t_ind(1):t_ind(2));
dt1 = yr_dt(t_ind(1):t_ind(2));

% Plot what this looks like
hFig = figure(1);
% set(hFig, 'Position', [100 100 920 500])
set(hFig, 'Position', [100 100 600 500])
hold on
plot(dt1, ld1, 'LineWidth', 2)
plot(dt1, pv1, 'LineWidth', 2)
ylabel("Power [kWac]")
legend("Load", "Solar PV", "location", "best", "interpreter", "latex")
% legend("Load Consumption", "Solar PV Generation", "location", "bestoutside")
set(gca, "FontSize", 28)
set(gca,'TickLabelInterpreter','latex');

%% 1st simple case: Fixed ESS size; fixed week, fixed initial costate; fixed intial SOC
% Get resulting curtailed PV, DG generation, final energy at optimal
% control given costate

% Costate
lambda_init = -80;

[u_opt, x, lambda] = mgpmpecm(length(ld1), pv1, ld1, x_max, x_min, LIB_INV_SIZE_KW, LIB_EFF_CHG, ...
    socdot, P_batt_range_socdot, SOC_range_socdot, ...
    dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
    h, lambda_init);

%WHY DOES U_OPT SOMETIMES EVOLVE SEPARATELY FROM SOC?? <-- due to NaNs

dg1 = max(ld1-pv1-u_opt,0);
pv_curt = max(pv1-ld1-u_opt,0);
x_f = x(end);
%%
hFig = figure();
set(hFig, 'Position', [100 100 1000 1000])
hold on
plot(dt1, ld1, 'LineWidth', 2)
plot(dt1, pv1, 'LineWidth', 2)
plot(dt1, u_opt, 'LineWidth', 3)
plot(dt1, dg1, 'LineWidth', 2)
ylabel("Power [kWAC]")
legend("Load", "PV", "Optimal LIB Dispatch", "Diesel Genset", "location", "eastoutside", 'interpreter', 'latex')
set(gca, "FontSize", 28)
set(gca,'TickLabelInterpreter','latex');

%% Plot Optimal Dispatch!!
hFig = figure();
set(hFig, 'Position', [100 100 2000 2200])
subplot(2,2,1)
hold on
plot(dt1, ld1, 'LineWidth', 1)
plot(dt1, pv1, 'LineWidth', 1)
plot(dt1, u_opt, 'LineWidth', 3)
plot(dt1, dg1, 'LineWidth', 1)
ylabel("Power [kWAC]")
legend("Load", "PV", "Optimal ESS Dispatch", "Diesel Genset", "location", "best")
set(gca, "FontSize", 20)
subplot(2,2,2)
plot(yr_dt(t_ind(1):t_ind(2)+1), x, 'LineWidth', 2)
ylim([0 1])
ylabel("SOC [-]")
set(gca, "FontSize", 20)
subplot(2,2,3)
plot(yr_dt(t_ind(1):t_ind(2)+1), lambda, 'LineWidth', 2)
ylabel("$\lambda$ [L]")
set(gca, "FontSize", 20)
subplot(2,2,4)
plot(dt1, cumsum(genset_model(dg1))*h, 'LineWidth', 2)
ylabel('Cumulative Fuel Consumption [L]')
set(gca, "FontSize", 20)
disp(sum(genset_model(dg1)*h))


%% 2nd simple case: Variation of initial costate
% Get resulting curtailed PV, DG generation, final energy at optimal
% control given costate

% Costate
lambda_inits = linspace(-150, -50, 41);

for i = 1:length(lambda_inits)
    [u_opt, x, lambda] = mgpmpecm(length(ld1), pv1, ld1, x_max, x_min, LIB_INV_SIZE_KW, LIB_EFF_CHG, ...
        socdot, P_batt_range_socdot, SOC_range_socdot, ...
        dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
        h, lambda_inits(i));
    u_opts(:,i) = u_opt;
    x_opts(:,i) = x;
    lambda_opts(:,i) = lambda;
    xfs(i) = x(end);
    dgi = max(ld1-pv1-u_opt,0);
    fuel_conss(i) = sum(genset_model(dgi))*h;
    disp(['Iteration complete: ', num2str(i)])
end
disp(['All ' num2str(i) ' iterations complete!'])

% Plot results of initial costate variation
% create a default color map ranging from red to blue
len = length(lambda_inits);
red = [1, 0, 0];
blue = [0, 1, 0];
colors_p = [linspace(red(1),blue(1),len)', linspace(red(2),blue(2),len)', linspace(red(3),blue(3),len)'];

hFig = figure(10);
set(hFig, 'Position', [100 100 800 700])
hold on
for i = 1:length(lambda_inits)
    plot(dt1, u_opts(:,i), 'Color', colors_p(i,:), 'LineWidth', 2)
end
% plot(u_opts, 'LineWidth', 2)
ylabel("Power [kWAC]")
title("Optimal LIB Dispatch vs. Initial Costate $\lambda_0$")
legendStrings = "\lambda_0 = " + string(lambda_inits);
legend(legendStrings, 'location', 'eastoutside')
set(gca, "FontSize", 20)

hFig = figure(11);
set(hFig, 'Position', [100 100 800 700])
hold on
for i = 1:length(lambda_inits)
    plot(yr_dt(t_ind(1):t_ind(2)+1), x_opts(:,i), 'Color', colors_p(i,:), 'LineWidth', 2)
end
ylabel("SOC [-]")
title("Optimal LIB SOC vs. Initial Costate $\lambda_0$")
% plot(x_opts, 'LineWidth', 2)
legendStrings = "\lambda_0 = " + string(lambda_inits);
legend(legendStrings, 'location', 'eastoutside')
set(gca, "FontSize", 20)

hFig = figure(12);
set(hFig, 'Position', [100 100 800 700])
hold on
for i = 1:length(lambda_inits)
    plot(yr_dt(t_ind(1):t_ind(2)+1), lambda_opts(:,i), 'Color', colors_p(i,:), 'LineWidth', 2)
end
% plot(lambda_opts, 'LineWidth', 2)
title("Evolution of Optimal $\lambda$ vs. Initial Costate $\lambda_0$")
ylabel("$\lambda$ [L]")
legendStrings = "\lambda_0 = " + string(lambda_inits);
legend(legendStrings, 'location', 'eastoutside')
set(gca, "FontSize", 20)

hFig = figure(13);
set(hFig, 'Position', [100 100 800 700])
hold on
for i = 1:length(lambda_inits)
    plot(lambda_inits(i), xfs(i), 'Color', colors_p(i,:), 'marker','^', 'MarkerSize',10, 'MarkerFaceColor', colors_p(i,:))
end
title("Final LIB SOC vs. Initial Costate $\lambda_0$")
xlabel("$\lambda$ [L]")
ylabel("SOC [-]")
% plot(lambda_opts, 'LineWidth', 2)
legendStrings = "\lambda_0 = " + string(lambda_inits);
legend(legendStrings, 'location', 'eastoutside')
set(gca, "FontSize", 20)

hFig = figure(14);
set(hFig, 'Position', [100 100 800 700])
hold on
for i = 1:length(lambda_inits)
    plot(lambda_inits(i), fuel_conss(i), 'Color', colors_p(i,:), 'marker','^', 'MarkerSize',10, 'MarkerFaceColor', colors_p(i,:))
end
title("Total Fuel Consumption vs. Initial Costate $\lambda_0$")
xlabel("$\lambda$ [L]")
ylabel("Diesel Consumption [L]")
% plot(lambda_opts, 'LineWidth', 2)
legendStrings = "\lambda_0 = " + string(lambda_inits);
legend(legendStrings, 'location', 'eastoutside')
set(gca, "FontSize", 20)

hFig = figure(15);
set(hFig, 'Position', [100 100 800 700])
hold on
for i = 1:length(lambda_inits)
    plot(xfs(i), fuel_conss(i), 'Color', colors_p(i,:), 'marker','^', 'MarkerSize',10, 'MarkerFaceColor', colors_p(i,:))
end
title("Total Fuel Consumption vs. Final LIB SOC")
xlabel("SOC [-]")
ylabel("Diesel Consumption [L]")
% plot(lambda_opts, 'LineWidth', 2)
legendStrings = "\lambda_0 = " + string(lambda_inits);
legend(legendStrings, 'location', 'eastoutside')
set(gca, "FontSize", 20)


%% SCRATCH

% yr_dt_sec = (YEAR_START:seconds(1):datetime(2021,12,31,23,59,59))';

% week_len = 4*24*7;
% day_len = 4*24;
% % week_start = 0
% week_start = 30 * week_len; % try 25?
% week_end = week_start + week_len-1;
% 
% day_start = week_start;
% day_end = day_start + day_len-1;

% % One week:
% ld_wk1 = ld(week_start:week_end);
% pv_wk1 = pv(week_start:week_end);
% wk1_dt = yr_dt(week_start:week_end);
% 
% % One day (START WITH THIS):
% ld1 = ld(day_start:day_end);
% pv1 = pv(day_start:day_end);
% d1_dt = yr_dt(day_start:day_end);

% % Convert to seconds
% ld_sec = repelem(ld1, 3600*h);
% pv_sec = repelem(pv1, 3600*h);
% d1_dt_sec = yr_dt_sec(30*24*7*3600+1:30*24*7*3600+(3600*24));

% %% 1st simple case: Fixed energy, power of ESS; fixed week, fixed costate; fixed starting energy
% % Initial value of costate
% lambda_init = 10;
% 
% [u_opt, x, lambda] = mgpmp(length(ld_sec), pv_sec, ld_sec, x_max, x_min, LIB_INV_SIZE_KW, LIB_EFF_CHG, ...
%     socdot, P_batt_range_socdot, SOC_range_socdot, ...
%     dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
%     h, lambda_init);
% 
% % Get resulting curtailed PV, DG generation, final energy at optimal
% % control given costate
% dg1 = max(ld_sec-pv_sec-u_opt,0);
% pv_curt = max(pv_sec-ld_sec-u_opt,0);
% x_f = x(end);
% 
% % Plot Optimal Dispatch!!
% hFig = figure(3);
% set(hFig, 'Position', [100 100 2000 2200])
% subplot(2,2,1)
% hold on
% plot(d1_dt_sec, ld_sec, 'LineWidth', 2)
% plot(d1_dt_sec, pv_sec, 'LineWidth', 2)
% plot(d1_dt_sec, u_opt, 'LineWidth', 2)
% plot(d1_dt_sec, dg1, 'LineWidth', 2)
% ylabel("Power [kWAC]")
% legend("Load", "PV", "Optimal ESS Dispatch", "Diesel Genset", "location", "best")
% set(gca, "FontSize", 20)
% subplot(2,2,2)
% plot(yr_dt_sec(30*24*7*3600+1:30*24*7*3600+(3600*24)+1), x, 'LineWidth', 2)
% ylim([0 1])
% ylabel("SOC [-]")
% set(gca, "FontSize", 20)
% subplot(2,2,3)
% plot(yr_dt_sec(30*24*7*3600+1:30*24*7*3600+(3600*24)+1), lambda, 'LineWidth', 2)
% ylabel("\lambda [L]")
% set(gca, "FontSize", 20)
% subplot(2,2,4)
% plot(d1_dt_sec, cumsum(genset_model(dg1))*h/3600, 'LineWidth', 2)
% ylabel('Cumulative Fuel Consumption [L]')
% set(gca, "FontSize", 20)

%% Function to calculate optimal control, final state, and state vector
% % TODO: DELETE THIS HERE, USE mgpmpecm.m INSTEAD
% function [u_opt, x, lambda] = mgpmp(outage_len, pv_outage, ld_outage, ...
%     x_max, x_min, inv_rating, inv_eff, ...
%     socdot, P_batt_range_socdot, SOC_range_socdot, ...
%     dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
%     h, lambda_init)
% % Function to compute optimal control strategy based on Pontryagin
% % Minimization Principle (PMP)
% % Inputs:
% % outage_len = length of outage period (e.g. one day) in 15 min. intervals
% % pv_outage = pv profile during outage period, as kW
% % ld_outage = load profile during outage period, as kW
% % x_max = upper bound of state
% % x_min = lower bound of state
% % inv_rating = kWAC rating of LIB inverter
% % inv_eff = AC/DC conversion efficiency of LIB inverter
% % socdot = map from (pack power, SOC) --> change in SOC
% % P_batt_range_socdot = corresponding pack power range for socdot map
% % SOC_range_socdot = corresponding SOC range for socdot map
% % dfdx = map from (pack power, SOC) --> dSOC_dot/dSOC
% % P_batt_range_dfdx = corresponding pack power range for dfdx map
% % SOC_range_dfdx = corresponding SOC range for dfdx map
% % lambda = initial costate
% % h = fraction of the hour
% % Outputs: 
% % u_opt = optimal control
% % x = state as a function of time from t=1:day_len+1
% % lambda = costate
% x = zeros(outage_len+1, 1);
% u_opt = zeros(outage_len, 1);
% lambda = zeros(outage_len+1, 1);
% 
% % TODO: Change this assumption
% x(1) = x_max;
% 
% lambda(1) = lambda_init;
% 
% % Set penalty (K) for violating constraints
% K = 200; % TODO: Sensitivity analysis on this
% 
% %     for t = 1:3600*5
%     for t = 1:outage_len
% %         disp('Time')
% %         disp(t)
%         % First, work in all AC power
%         % Determine control space at time t based on PV, load power
%         % assuming that ESS can only charge off of PV, discharge to load
%         % But bounded by inverter rating
%         u_min_t = -min(pv_outage(t), inv_rating);
%         u_max_t = min(ld_outage(t), inv_rating);
% %         u_min_t = -pv_outage(t);
% %         u_max_t = ld_outage(t);
%         u_range_t = sort([linspace(u_min_t,u_max_t) 0]); % include 0 power as an option
% %         disp(u_range_t)
% 
%         % Calculate optimal control with control space at time t
%         H_t = zeros(length(u_range_t), 1);
%         wE_t = zeros(length(u_range_t), 1);
% 
%         for v = 1:length(u_range_t)
%             % Put rules (~ environment) for each dispatch here
%             % TODO: Break this out into its own function (or file/method??)
%             % Initialize PV, load, and ESS energy resources; DG power, penalty
%             pv_t = pv_outage(t);
%             l_t = ld_outage(t);
%             x_t = x(t);
%             dg_t = 0;
% 
%             u_t = u_range_t(v);
% 
%             % First, charge PV off of ESS if u is negative
%             pv_t = pv_t - (-min(0,u_t));
% 
%             % Dispatch portion of load from ESS if u is positive
%             l_t = l_t - (max(0,u_t));
% 
%             % Balance the load and PV, with excess PV curtailed or excess load
%             % supplied by DG
%             if pv_t > l_t
%                 % Subtract off load from PV and curtail the rest
%                 % DG is not used here
%             elseif pv_t < l_t
%                 % Subtract off PV from load and supply remainder load with DG
%                 l_t = l_t - pv_t;
%                 dg_t = l_t;
%             else
%                 % Load and PV are in perfect balance! How??
%                 % DG is not used here
%             end
% %             disp(dg_t)
%             % Calculate penalty contribution
%             % First, convert AC power to DC power
%             if u_t >= 0
%                 u_DC_t = u_t/inv_eff;
%             else
%                 u_DC_t = u_t * inv_eff;
%             end
%             
%             % Use socdot map to find SOC at time step t+1
%             socdot_t = interp2(P_batt_range_socdot, SOC_range_socdot, socdot, u_DC_t, x_t, 'spline');
%             
%             x_t1 = x_t + socdot_t*h;
% %             disp(x_t1)
% %             x_t1 = x_t + socdot_t;
%             if x_t1 > x_max
%                 w = K;
%             elseif x_t1 < x_min
%                 w = -K;
%             elseif x_min <= x_t1 && x_t1 <= x_max
%                 w = 0;
% %             else
% %                 w = K;
%             end
% %             disp(w)
%             wE_t(v) = w;
% %             disp(w)
%             % Calculate fuel consumpt ion from function
%             fuel_cons_L = genset_model(dg_t); 
% %             fuel_cons_L = genset_model(dg_t) * h / 3600; 
% 
%             % Put resulting Hamiltonian here
% %             H_t(v) = fuel_cons_L + (lambda(t) + w) * u_t;
%             H_t(v) = fuel_cons_L + (lambda(t) + w) * socdot_t;
% %             disp(H_t(v))
%         end
% %         % Plot
% %         figure(2)
% %         subplot(2,1,1)
% %         plot(u_range_t,H_t)
% %         ylabel("Hamiltonian value")
% %         xlabel("Power, kW")
% %         subplot(2,1,2)
% %         plot(u_range_t, wE_t, 'o')
% %         ylabel("Penalty value")
% %         xlabel("Power, kW")
%         
%         % Get index of H_t where it is minimized;
%         [~, min_ind] = min(H_t);
% %         disp(min(H_t))
% %         disp(min_ind)
%         % Find the corresponding LIB AC power value as the optimal control
%         u_opt_t = u_range_t(min_ind);
%         
%         % Store optimal LIB AC power control
%         u_opt(t) = u_opt_t;
% 
%         % First, convert AC power to DC power
%         if u_opt_t >= 0
%             u_DC_opt_t = u_opt_t/inv_eff;
%         else
%             u_DC_opt_t = u_opt_t * inv_eff;
%         end
% 
%         % Use socdot map to find SOC at time step t+1
%         socdot_t = interp2(P_batt_range_socdot, SOC_range_socdot, socdot, u_DC_opt_t, x(t), 'spline');
%         
%         % Update the SOC for the next time step
% %         x(t+1) = x(t) + socdot_t;
%         x(t+1) = x(t) + socdot_t*h;
% %         disp(x(t+1))
%         % Update the costate
%         % 1) Calculate penalty
%         if x(t) > x_max
%             w = K;
%         elseif x(t) < x_min
%             w = -K;
%         elseif x_min <= x(t) && x(t) <= x_max
%             w = 0;
% %         else
% %             w = K;
%         end
%         
%         % 2) Compute dfdx
%         dfdx_t = -interp2(P_batt_range_dfdx, SOC_range_dfdx, dfdx, u_DC_opt_t, x(t), 'spline');
%         
%         % 3) Compute lambda_dot
%         lambda_dot = -(lambda(t) + w)*dfdx_t;
% 
%         % 4) Evolve costate
%         lambda(t+1) = lambda(t) + lambda_dot*h;
%     end
% 
%     
%     % function end
% end
