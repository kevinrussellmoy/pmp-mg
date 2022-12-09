%% File to generate diesel genset data
% Source of data: https://akenergygateway.alaska.edu/media/AMP/AMP%20General/Genator%20Fuel%20Consumpion%20Under%20Dynamic%20Loading%20201710.pdf
% Kevin Moy
% Sep 8 2021
clearvars
close all
clc
%% Load data obtained via GRABIT function from Fig. 12
load('gen_fuel_curve.mat');

load_kw = gen_fuel_curve(:,1); % in kW
fuel_eff = gen_fuel_curve(:,2); % in kWh/L

fuel_cons = load_kw ./ fuel_eff; % in L/h

% try removing 1st data pt as an outlier

load_kw2 = load_kw(2:end);
fuel_cons2 = fuel_cons(2:end);

%% Fit to 2nd order polynomial

B2 = polyfit(load_kw2, fuel_cons2, 2);
% RMSE 0.3242 (rsquare 0.9999)
% 1.2349e-4 load^2 + 0.1982 load + 16.3602 [L/hr]
plot(load_kw2, fuel_cons2)

hFig = figure(1);
plot(load_kw2, fuel_cons2, 'LineWidth', 2)
xlabel('Genset Output [kWac]', 'Interpreter','latex')
ylabel('Fuel consumption [L/hr]', 'Interpreter','latex')
% title('$\partial \dot{SOC} / \partial SOC$', 'Interpreter','latex')
ax = gca;
set(ax,'FontSize',28);
box on