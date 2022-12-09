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
% set(hFig, 'Position', [100 100 1000 500])
set(hFig, 'Position', [100 100 600 600])
hold on
plot(dt1, ld1, 'LineWidth', 2)
plot(dt1, pv1, 'LineWidth', 2)
plot(dt1, dg1, 'LineWidth', 2)
plot(dt1, u_opt, 'LineWidth', 5, 'color', rgb('green'))
ylabel("Power [kWAC]")
% legend("Load", "PV", "Diesel Genset", "Optimal LIB Dispatch", "location", "eastoutside", 'interpreter', 'latex')
legend("Load", "PV", "Diesel Genset", "Optimal LIB Dispatch", "location", "southoutside", 'interpreter', 'latex')
set(gca, "FontSize", 28)
set(gca,'TickLabelInterpreter','latex');

%% Plot Optimal Dispatch!!
hFig = figure();
set(hFig, 'Position', [100 100 2000 2200])
subplot(2,2,1)
hold on
plot(dt1, ld1, 'LineWidth', 1)
plot(dt1, pv1, 'LineWidth', 1)
plot(dt1, u_opt, 'LineWidth', 5, 'color', rgb('green'))
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
lambda_inits = linspace(-150, -50, 21);

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

