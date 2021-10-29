%% PMP Sizing for Microgrid Energy Storage
% Kevin Moy, 8/15(?)/21
clearvars
close all

%% Global vartiables(sizes):
PV_ARRAY_SIZE_KW = 660;   % kWAC rating of the PV array
DIESEL_GEN_SIZE_KW = 925.0;   % kWAC rating of the diesel generator
% % Diesel fuel consumption coefficients from https://ieeexplore.ieee.org/document/8494571
% DIESEL_FUEL_CONS_A =   0.246; % Liters per kWh
% DIESEL_FUEL_CONS_B =   0.08415; % Liters per kW (rating)

STORAGE_DURATION = 8; % Hours of storage duration at maximum power

ESS_EFF_DISCHG = 0.95; % Efficiency of discharging ESS
ESS_EFF_CHG = 0.95; % Efficiency of charging ESS

% Read in CSV data sources:
bldg = readmatrix('bldg_load.csv');
pv_r = readmatrix('pv_gen.csv');
resi = readmatrix('resi_load.csv');

%% Diesel fuel consumption table from: https://ieeexplore.ieee.org/abstract/document/4276041

loadpct = [0.2; 0.25; 0.3; 0.4; 0.5; 0.6; 0.7; 0.75; 0.8; 0.9; 1.0];
fuelrate = [61.4; 71.9; 82.3; 102.8; 123.0; 144.0; 165.2; 175.9; 186.6; 209.2; 233.0];
fueleff = [3.01; 3.22; 3.37; 3.60; 3.76; 3.85; 3.92; 3.92; 3.97; 3.98; 3.97];

hFig = figure(1);
set(hFig, 'Position', [100 100 600 500])
hold on
yyaxis left
plot(loadpct, fuelrate, 'o', 'LineWidth', 2, 'MarkerSize',10)
xlabel('Percent of Rated Power')
ylabel('Fuel Rate (L/hr)')
ylim([0 240])
yyaxis right
plot(loadpct, fueleff, 'x', 'LineWidth', 2, 'MarkerSize',10)
ylabel('Fuel Efficiency (kWh/L)')
ylim([0 4])
set(gca, "FontSize", 20)
xlim([0 1])

%% Datetime vector for plotting
YEAR_START = datetime(2021,1,1,0,0,0);
YEAR_END = datetime(2021,12,31,23,45,0);
yr_dt = (YEAR_START:minutes(15):YEAR_END)';

%% Obtain aggregate load:
load = bldg(:,2)*0.25 + bldg(:,3)*0.25 + resi(:,2);
pv = pv_r(:,2) * PV_ARRAY_SIZE_KW / max(pv_r(:,2));

%% Plot load, PV data with datetime x-axis
hFig = figure(1);
set(hFig, 'Position', [100 100 600 500])
hold on
plot(yr_dt, load)
plot(yr_dt, pv)

% Plot residual load per day
hFig = figure(2);
set(hFig, 'Position', [100 100 600 500])
hold on
plot(yr_dt, load-pv)

%% 
% ''' ~.~.~.~ optimization time ~.~.~.~ '''
% randomly (or not so randomly) select 7-day intervals to optimize the dispatch
% first do for a set ESS size (500 kW, 950 kWh as in BLR Microgrid)
% then make the ESS size a part of the function!
% Constrain storage size to [min, max] and similarly power
% Then find the optimum, and then find the closest "round" value and present those power flows
% then add in degradation penalty


week_len = 4*24*7;
day_len = 4*24;
% week_start = 0
week_start = 30 * week_len;
week_end = week_start + week_len-1;

day_start = week_start;
day_end = day_start + day_len-1;

% One week:
ld_wk1 = load(week_start:week_end);
pv_wk1 = pv(week_start:week_end);
wk1_dt = yr_dt(week_start:week_end);

% One day (START WITH THIS):
ld1 = load(day_start:day_end);
pv1 = pv(day_start:day_end);
d1_dt = yr_dt(day_start:day_end);

% Fraction of the hour
h = 15/60;

% Plot what this looks like
figure(2)
hold on
plot(wk1_dt, ld_wk1)
plot(wk1_dt, pv_wk1)
figure(3)
hold on
plot(d1_dt, ld1)
plot(d1_dt, pv1)

%% 1st simple case: Fixed energy, power of ESS; fixed week, fixed costate; fixed starting energy

% Rated Energy, Rated Power, Energy limits of the ESS
E_nom = 1000; % kWh
P_nom = 500; % kW
% TODO consider different min/max energy values?
E_min = 0;
E_max = E_nom;

% Set penalty (K) for violating constraints
K = 100; % TODO: Sensitivity analysis on this

% Store relvant quantities
% TODO: Squash some of these into arrays, e.g. PV power flow into a 3xN
% matrix
E = zeros(day_len+1, 1);
H = zeros(day_len, 1);
u_opt = zeros(day_len, 1);
% TODO: Calculate these and store
pv_ess = zeros(day_len, 1);
pv_l = zeros(day_len, 1);
pv_curt = zeros(day_len, 1);
ess_l = zeros(day_len, 1);
% l_curt = zeros(day_len, 1); % For now, no load curtailment; so that all
% supplied load must be met by PV, ESS, and DG
dg = zeros(day_len, 1);

% Costate
lambda = 8;

% % Sweep of control space, where negative is charge, positive is discharge
% u_range = -P_nom:10:P_nom;

% Set initial stored energy:
E(1) = E_nom; % TODO: Relax this assumption later

[u_opt, E] = mgpmp(day_len, pv1, ld1, E_nom, 0, h, lambda);

% Get resulting curtailed PV, DG generation, final energy at optimal
% control given costate
dg1 = max(ld1-pv1-u_opt,0);
pv_curt = max(pv1-ld1-u_opt,0);
E_f = E(day_len+1);

% Plot Optimal Dispatch!!

hFig = figure(3);
set(hFig, 'Position', [100 100 800 1000])
subplot(2,1,1)
hold on
plot(d1_dt, ld1, 'LineWidth', 2)
plot(d1_dt, pv1, 'LineWidth', 2)
plot(d1_dt, u_opt, 'LineWidth', 2)
plot(d1_dt, dg1, 'LineWidth', 2)
ylabel("Power, kW")
legend("Load", "PV", "Optimal ESS Dispatch", "Diesel Genset", "location", "best")
set(gca, "FontSize", 20)
subplot(2,1,2)
plot(yr_dt(day_start:day_end+1), E, 'LineWidth', 2)
ylim([E_min, 1.1*E_max])
ylabel("Stored ESS Energy, kWh")
set(gca, "FontSize", 20)

% TODO: Compare this to Gurobi??

%% 2nd simple case: 1st simple case but sweep on costate lambda for sensitivity analysis

lambdas = linspace(0,15,16);
E_finals = zeros(length(lambdas),1);
tic
for i = 1:length(lambdas)
    lambda = lambdas(i);
    [~, E] = mgpmp(day_len, pv1, ld1, E_nom, 0, h, lambda);
    E_finals(i) = E(end);
end
toc

hFig = figure(10);
set(hFig, 'Position', [100 100 600 500])
plot(lambdas, E_finals, 'o', 'LineWidth', 2, 'MarkerSize',10)
xlabel('Costate \lambda_0')
ylabel('Final state E(t_f)')
set(gca, "FontSize", 20)


%% 2nd simple case starting at 200 random 15-minute intervals each year
% TODO: do for each 15-minute period in the year, but it takes too long
% ~1.25 seconds * (35040-(4*24)) = big number
% E_finals_year = zeros(length(load)-day_len-1, length(lambdas));

% initialize random number generator
rng(0,'twister');

% Generate 200 random outage starting periods and sort (so datetimes are in
% order)
num_outages = 200;
r = randi([1 length(load)-day_len-1],num_outages,1);
r_sort = sort(r);

% Store final state for each random starting period
E_finals_rand = zeros(length(r_sort), length(lambdas));

% Get datetimes for each starting period
yr_dt_rand = yr_dt(r_sort);


% for i = 1:length(load)-day_len-1
tic
for i = 1:num_outages
    ld_day = load(i:i+day_len);
    pv_day = pv(i:i+day_len);
    for j = 1:length(lambdas)
        lambda = lambdas(j);
        [~, E] = mgpmp(day_len, pv_day, ld_day, E_nom, 0, h, lambda);
        E_finals_rand(i,j) = E(end);
    end
end
toc
%% plotting above
hFig = figure(10);
set(hFig, 'Position', [100 100 600 500])
surf(lambdas, yr_dt_rand, E_finals_rand)
ylabel('Starting time of outage')
xlabel('Costate \lambda_0')
zlabel('Final state E(t_f), kWh')
title({'1 day of outage from', '200 random starting 15-minute periods'})
set(gca, "FontSize", 20)
ax = gca;
set(ax,'FontSize',28);
cb = colorbar('Location','eastoutside');
ax.Position = ax.Position - [0 0 .1 .1];
cb.Position = cb.Position + [.1 0 0 0];

% color scheme
% c = turbo(200);
c = parula(length(load));

hFig = figure(11);
set(hFig, 'Position', [100 100 600 500])
hold on
for i = 1:num_outages
    c_i = c(r_sort(i),:);
    plot(lambdas, E_finals_rand(i,:), 'Color', c_i)
end
colorbar
ylabel('Starting time of outage')
xlabel('Costate \lambda_0')
set(gca, "FontSize", 20)


%% 3rd simple case: 1st simple case but with final state determination via bisectional search for costate lambda






%% Function to calculate optimal control, final state, and state vector

function [u_opt, E] = mgpmp(outage_len, pv_outage, ld_outage, E_max, E_min, h, lambda)
% Function to compute optimal control strategy based on Pontryagin
% Minimization Principle (PMP)
% Inputs:
% outage_len = length of outage period (e.g. one day) in 15 min. intervals
% pv_outage = pv profile during outage period
% ld_outage = load profile during outage period
% E_max = upper bound of state
% E_min = lower bound of state
% lambda = costate
% h = fraction of the hour
% Outputs: 
% u_opt = optimal control
% E = state as a function of time from t=1:day_len+1
E = zeros(outage_len+1, 1);
u_opt = zeros(outage_len, 1);
E(1) = E_max;

% Diesel fuel consumption table from: https://ieeexplore.ieee.org/abstract/document/4276041
DIESEL_GEN_SIZE_KW = 925.0;   % kWAC rating of the diesel generator
loadpct = [0.2; 0.25; 0.3; 0.4; 0.5; 0.6; 0.7; 0.75; 0.8; 0.9; 1.0];
fuelrate = [61.4; 71.9; 82.3; 102.8; 123.0; 144.0; 165.2; 175.9; 186.6; 209.2; 233.0];

% Set penalty (K) for violating constraints
K = 100; % TODO: Sensitivity analysis on this

    for t = 1:outage_len
        % Determine control space at time t based on PV, load power
        % and assuming that ESS can only charge off of PV, discharge to load
        u_min_t = -pv_outage(t);
        u_max_t = ld_outage(t);
        u_range_t = [linspace(u_min_t,u_max_t) 0]; 
        % TODO: Refine this range so that it is whole numbers, 
        % or at least always includes zero

        % Calculate optimal control with control space at time t
        H_t = zeros(length(u_range_t), 1);
        wE_t = zeros(length(u_range_t), 1);

        for v = 1:length(u_range_t)
            % Put rules (~ environment) for each dispatch here
            % TODO: Break this out into its own function (or file/method??)
            % Initialize PV, load, and ESS energy resources; DG power, penalty
            pv_t = pv_outage(t);
            l_t = ld_outage(t);
            E_t = E(t);
            dg_t = 0;
            w_E = 0;

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
            E_t = E_t - u_t*h;
            if E_t > E_max
                w_E = -K;
            elseif E_t < E_min
                w_E = K;
            else
                w_E = 0;
            end
            wE_t(v) = w_E;            

            % Put resulting Hamiltonian here
            H_t(v) = interp1(loadpct, fuelrate, dg_t/DIESEL_GEN_SIZE_KW, 'linear','extrap') * h * dg_t + (lambda + w_E) * u_t;
        end
    %     figure(2)
    %     subplot(2,1,1)
    %     plot(u_range_t,H_t)
    %     ylabel("Hamiltonian value")
    %     xlabel("Power, kW")
    %     subplot(2,1,2)
    %     plot(u_range_t, wE_t, 'o')
    %     ylabel("Penalty value")
    %     xlabel("Power, kW")
        % Get index of H_t where it is minimized;
        [~, min_ind] = min(H_t);
        % Find the corresponding u_vals as the optimal control
        u_opt_t = u_range_t(min_ind);
        u_opt(t) = u_opt_t;

        % Update the stored energy for the next time step
        E(t+1) = E(t) - h * u_opt_t;

    end

end


%% SCRATCH

% for t = 1:day_len

% % Function to compute optimal control strategy based on Pontryagin
% % Minimization Principle (PMP)
% % Inputs:
% % outage_len = length of outage period (e.g. one day) in 15 min. intervals
% % pv_outage = pv profile during outage period
% % ld_outage = load profile during outage period
% % Outputs: 
% % u_opt = optimal control
% % E = state as a function of time from t=1:day_len+1
% 
% for t = 1:day_len
%     % Determine control space at time t based on PV, load power
%     % and assuming that ESS can only charge off of PV, discharge to load
%     u_min_t = -pv1(t);
%     u_max_t = ld1(t);
%     u_range_t = [linspace(u_min_t,u_max_t) 0]; 
%     % TODO: Refine this range so that it is whole numbers, 
%     % or at least always includes zero
%     
%     % Calculate optimal control with control space at time t
%     H_t = zeros(length(u_range_t), 1);
%     wE_t = zeros(length(u_range_t), 1);
%     
%     for v = 1:length(u_range_t)
%         % Put rules (~ environment) for each dispatch here
%         % TODO: Break this out into its own function (or file/method??)
%         % Initialize PV, load, and ESS energy resources; DG power, penalty
%         pv_t = pv1(t);
%         l_t = ld1(t);
%         E_t = E(t);
%         dg_t = 0;
%         w_E = 0;
%         
%         u_t = u_range_t(v);
%         
%         % First, charge PV off of ESS if u is negative
%         pv_t = pv_t - (-min(0,u_t));
%         
%         % Dispatch portion of load from ESS if u is positive
%         l_t = l_t - (max(0,u_t));
%         
%         % Balance the load and PV, with excess PV curtailed or excess load
%         % supplied by DG
%         if pv_t > l_t
%             % Subtract off load from PV and curtail the rest
%             % DG is not used here
%         elseif pv_t < l_t
%             % Subtract off PV from load and supply remainder load with DG
%             l_t = l_t - pv_t;
%             dg_t = l_t;
%         else
%             % Load and PV are in perfect balance! How??
%             % DG is not used here
%         end
%                 
%         % Calculate penalty contribution
%         E_t = E_t - u_t*h;
%         if E_t > E_max
%             w_E = -K;
%         elseif E_t < E_min
%             w_E = K;
%         else
%             w_E = 0;
%         end
%         wE_t(v) = w_E;            
%         
%         % Put resulting Hamiltonian here
%         H_t(v) = interp1(loadpct, fuelrate, dg_t/DIESEL_GEN_SIZE_KW, 'linear','extrap') * h * dg_t + (lambda + w_E) * u_t;
%     end
% %     figure(2)
% %     subplot(2,1,1)
% %     plot(u_range_t,H_t)
% %     ylabel("Hamiltonian value")
% %     xlabel("Power, kW")
% %     subplot(2,1,2)
% %     plot(u_range_t, wE_t, 'o')
% %     ylabel("Penalty value")
% %     xlabel("Power, kW")
%     % Get index of H_t where it is minimized;
%     [minH, min_ind] = min(H_t);
%     % Find the corresponding u_vals as the optimal control
%     u_opt_t = u_range_t(min_ind);
%     u_opt(t) = u_opt_t;
%         
%     % Update the stored energy for the next time step
%     E(t+1) = E(t) - h * u_opt_t;
%     
% end



