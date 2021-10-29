%% Bisection search for optimal initial costate
clearvars
close all
clc
set(0,'defaultTextInterpreter','latex');
%% Load in microgrid configuration
mg_config

%% Define initial, final time

t_ini_str = '09-Oct-2019 00:00:00';

t_fin_str = '12-Oct-2019 17:45:00';

t_vec = [YEAR_START; datetime([t_ini_str; t_fin_str])];

t_ind = ((minutes(t_vec(2:3) - t_vec(1)))/15) + 1;

%% Select load and PV data

ld1 = ld(t_ind(1):t_ind(2));
pv1 = pv(t_ind(1):t_ind(2));
dt1 = yr_dt(t_ind(1):t_ind(2));

% Create variable for outage length (v useful)
outage_len = length(ld1);


%% target final SOC

SOC_f_target = 0.5;

% Bounds for search on initial value of costate
lambda_max = 0;
lambda_min = -200;

[lambda_inits, SOC_fins] = init_costate_search(length(ld1), pv1, ld1, x_max, x_min, ...
                BESS_INV_SIZE_KW, BESS_EFF_CHG, ...
                socdot, P_batt_range_socdot, SOC_range_socdot, ...
                dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
                h, lambda_min, lambda_max, SOC_f_target);


%%
% deltaSOC_target = 0;
% deltaSOC_tol = 0.01; % tolerance for accepting solution
% 
% s_min = -150;
% s_max = 0;
% s_init = -100;
% MaxNumIter = 10;
% numIter = 0;
% 
% tic
% while numIter < MaxNumIter
% numIter = numIter+1;
%     
%     s_init = (s_min + s_max)/2;
%     
%     [~, SOC, ~] = mgpmpecm(length(ld1), pv1, ld1, x_max, x_min, BESS_INV_SIZE_KW, BESS_EFF_CHG, ...
%     socdot, P_batt_range_socdot, SOC_range_socdot, ...
%     dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
%     h, s_init);
%     
%     SOC_end = SOC(end);
%     deltaSOC = SOC_end - SOC_f_target;
%     
%     s_init_vec(numIter) = s_init;
%     deltaSOC_vec(numIter) = deltaSOC;
%     SOC_end_vec(numIter) = SOC_end;
%     
%     disp(abs(deltaSOC - deltaSOC_target))
%     if abs(deltaSOC - deltaSOC_target) < deltaSOC_tol
%         break  % exit "while" loop
%     elseif deltaSOC > deltaSOC_target
%         s_min = s_init;
%     elseif deltaSOC < deltaSOC_target
%         s_max = s_init;
%     end
%     
% end
% toc

%% Plot!

hFig = figure(1);
set(hFig, 'Position', [100 100 600 2000])
subplot(3,1,1)
plot(lambda_inits, '.k', 'MarkerSize', 25)
xlim([1 length(SOC_fins)])
ylabel("Initial Costate $\lambda$ [L/hr]")
set(gca, "FontSize", 28)
% 
% subplot(3,1,2)
% plot(deltaSOC_vec, '.k', 'MarkerSize', 25)
% xlim([1 length(s_init_vec)])
% ylabel("SOC delta [-]")
% set(gca, "FontSize", 28)

subplot(3,1,3)
plot(SOC_fins, '.k', 'MarkerSize', 25)
xlim([1 length(SOC_fins)])
xlabel("Iteration")
ylabel("Ending SOC [-]")
set(gca, "FontSize", 28)

hFig = figure(2);
set(hFig, 'Position', [100 100 600 500])
plot(lambda_inits, SOC_fins, 'Marker', '.', 'MarkerSize', 25, 'LineWidth', 2) 
xlabel("Initial Costate $\lambda$ [L/hr]")
ylabel("Ending SOC [-]")
set(gca, "FontSize", 28)




