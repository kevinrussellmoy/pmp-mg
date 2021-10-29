%% View open-circuit voltage, internal resistance as previously identified by Edoardo
% See: https://doi.org/10.1016/j.apenergy.2021.116473
% Kevin Moy
% 8/23/2021
clearvars
close all
%% NMC, all cells

NMC_num_cells = 6;

load('/Users/kmoy14/Documents/Research/Model parameter identification/Identification_results_chemistry_comparison/R0_iden_NMC25degC.mat')
load('/Users/kmoy14/Documents/Research/Model parameter identification/Identification_results_chemistry_comparison/Vocv_iden_NMC25degC.mat')


% Concatenate SOC (R0), R0, SOC (OCV), OCV
SOC_R0_NMC = {SOC_points_NMC25_R0_M1; SOC_points_NMC25_R0_M2; ...
    SOC_points_NMC25_R0_M3; SOC_points_NMC25_R0_M4; ...
    SOC_points_NMC25_R0_M5; SOC_points_NMC25_R0_M6};

R0_NMC = {R0_NMC25_M1; R0_NMC25_M2; R0_NMC25_M3; R0_NMC25_M4; ...
    R0_NMC25_M5; R0_NMC25_M6};

SOC_OCV_NMC = {SOC_points_NMC25_Vocv_M1; SOC_points_NMC25_Vocv_M2; ...
    SOC_points_NMC25_Vocv_M3; SOC_points_NMC25_Vocv_M4; ...
    SOC_points_NMC25_Vocv_M5; SOC_points_NMC25_Vocv_M6};

OCV_NMC = {V_points_NMC25_M1; V_points_NMC25_M2; V_points_NMC25_M3;...
    V_points_NMC25_M4; V_points_NMC25_M5; V_points_NMC25_M6};

hFig = figure(1);
set(hFig, 'Position', [100 100 800 1000])
subplot(2,1,1)
hold on
for i = 1:NMC_num_cells
    plot(SOC_R0_NMC{i}, R0_NMC{i}, 'LineWidth', 2)
end
legendCell = cellstr(num2str((1:NMC_num_cells)', 'Cell %-d'));
xlabel('SOC, [-]')
ylabel('R_0, [Ohm]')
title('Internal Resistance vs. SOC')
legend(legendCell, 'Location', 'eastoutside')

set(gca, "FontSize", 20)

subplot(2,1,2)
hold on
for i = 1:NMC_num_cells
    plot(SOC_OCV_NMC{i}, OCV_NMC{i}, 'LineWidth', 2)
end
legendCell = cellstr(num2str((1:NMC_num_cells)', 'Cell %-d'));
xlabel('SOC, [-]')
ylabel('V_{OC}, [V]')
title('Open Circuit Voltage vs. SOC')
legend(legendCell, 'Location', 'eastoutside')
sgtitle('NMC', 'FontSize', 22)
set(gca, "FontSize", 20)

%% LFP, all cells

LFP_num_cells = 6;

load('/Users/kmoy14/Documents/Research/Model parameter identification/Identification_results_chemistry_comparison/R0_iden_LFP25degC.mat')
load('/Users/kmoy14/Documents/Research/Model parameter identification/Identification_results_chemistry_comparison/Vocv_iden_LFP25degC.mat')


% Concatenate SOC (R0), R0, SOC (OCV), OCV
SOC_R0_LFP = {SOC_points_LFP25_R0_M1; SOC_points_LFP25_R0_M2; ...
    SOC_points_LFP25_R0_M3; SOC_points_LFP25_R0_M4; ...
    SOC_points_LFP25_R0_M5; SOC_points_LFP25_R0_M6};

R0_LFP = {R0_LFP25_M1; R0_LFP25_M2; R0_LFP25_M3; R0_LFP25_M4; ...
    R0_LFP25_M5; R0_LFP25_M6};

SOC_OCV_LFP = {SOC_points_LFP25_Vocv_M1; SOC_points_LFP25_Vocv_M2; ...
    SOC_points_LFP25_Vocv_M3; SOC_points_LFP25_Vocv_M4; ...
    SOC_points_LFP25_Vocv_M5; SOC_points_LFP25_Vocv_M6};

OCV_LFP = {V_points_LFP25_M1; V_points_LFP25_M2; V_points_LFP25_M3;...
    V_points_LFP25_M4; V_points_LFP25_M5; V_points_LFP25_M6};

hFig = figure(1);
set(hFig, 'Position', [100 100 800 1000])
subplot(2,1,1)
hold on
for i = 1:LFP_num_cells
    plot(SOC_R0_LFP{i}, R0_LFP{i}, 'LineWidth', 2)
end
legendCell = cellstr(num2str((1:LFP_num_cells)', 'Cell %-d'));
xlabel('SOC, [-]')
ylabel('R_0, [Ohm]')
title('Internal Resistance vs. SOC')
legend(legendCell, 'Location', 'eastoutside')

set(gca, "FontSize", 20)

subplot(2,1,2)
hold on
for i = 1:LFP_num_cells
    plot(SOC_OCV_LFP{i}, OCV_LFP{i}, 'LineWidth', 2)
end
legendCell = cellstr(num2str((1:LFP_num_cells)', 'Cell %-d'));
xlabel('SOC, [-]')
ylabel('V_{OC}, [V]')
title('Open Circuit Voltage vs. SOC')
legend(legendCell, 'Location', 'eastoutside')
sgtitle('LFP', 'FontSize', 22)
set(gca, "FontSize", 20)


