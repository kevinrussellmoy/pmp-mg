%% Calculate function approximations for NMC open-circuit voltage, internal resistance as previously identified by Edoardo
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
box on
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
box on
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

%% Clear unneeded variables
clearvars -except SOC_OCV_NMC OCV_NMC SOC_R0_NMC R0_NMC NMC_num_cells

%% Extract V_OC curves

% Keep only cells 1, 3, 6
keepcells = [1; 3; 6];

SOC_min = max([min(SOC_OCV_NMC{1}), min(SOC_OCV_NMC{3}), min(SOC_OCV_NMC{6})]);

% Interpolate so they're all at the same length
SOC_new = flip(linspace(SOC_min,1,30));

for i = 1:length(keepcells)
    SOC_OCV = SOC_OCV_NMC{keepcells(i)};
    OCV = OCV_NMC{keepcells(i)};
    OCV_new(:,i) = interp1(SOC_OCV, OCV, SOC_new, 'linear','extrap');
end

% hFig = figure(2);
% set(hFig, 'Position', [100 100 600 500])
% hold on
% xlabel('SOC, [-]')
% ylabel('V_{OC}, [V]')
% title('OCV vs. SOC')
% plot(SOC_new, OCV_new, 'LineWidth', 2)
% legend('Cell 1', 'Cell 3', 'Cell 6')
% box on
% set(gca, "FontSize", 24)


% Average to get the new OCV(SOC) curve
OCV_avg = mean(OCV_new,2);

% Two splines: first from SOC_new(1:4) and second from SOC_new(4:end)
SOC_1 = SOC_new(27:end);
OCV_1 = OCV_avg(27:end);
B1 = polyfit(SOC_1, OCV_1, 1);
F1 = polyval(B1, SOC_1);
RMSE1 = sqrt(mean((OCV_1' - F1).^2));
disp(RMSE1)
% Results:
% OCV_hat_1 = 2.6071 [V] * (SOC) + 2.9333 [V]
% RMSE: 0.0077
% for SOC in [0.0949, 0.1885]

SOC_2 = SOC_new(1:27);
OCV_2 = OCV_avg(1:27);
B2 = polyfit(SOC_2, OCV_2, 1);
F2 = polyval(B2, SOC_2);
RMSE2 = sqrt(mean((OCV_2' - F2).^2));
disp(RMSE2)
% Results:
% OCV_hat_2 = 0.9374 [V] * (SOC) + 3.2436 [V]
% RMSE: 0.0111
% for SOC in (0.1885, 1]

% Plot
hFig = figure(3);
set(hFig, 'Position', [100 100 600 500])
box on
hold on
plot(SOC_new, OCV_new, '.', 'MarkerSize', 15)
plot(SOC_1, F1, 'k', 'LineWidth', 2.5)
plot(SOC_2, F2, 'k', 'LineWidth', 2.5)
xlabel('SOC [-]')
ylabel('$V_{oc}$ [V]')
% title('Open-Circuit Voltage vs. SOC')
legend({'Cell 1', 'Cell 3', 'Cell 6', ['Averaged Piecewise' newline 'Linear Fit']}, 'Location', 'southeast', 'interpreter', 'latex')
% legend('','Spline 1','Spline 2', 'Location', 'southeast')
box on
set(gca, "FontSize", 28)

%% Extract R_0 curves

% Keep only cells 1, 3, 6
keepcells = [1; 3; 6];

SOC_min = max([min(SOC_R0_NMC{1}), min(SOC_R0_NMC{3}), min(SOC_R0_NMC{6})]);

% Interpolate so they're all at the same length
SOC_new = flip(linspace(SOC_min,1,30));

for i = 1:length(keepcells)
    SOC_R0 = SOC_R0_NMC{keepcells(i)};
    R0 = R0_NMC{keepcells(i)};
    R0_new(:,i) = interp1(SOC_R0, R0, SOC_new, 'linear','extrap');
end

% hFig = figure(4);
% set(hFig, 'Position', [100 100 600 500])
% hold on
% xlabel('SOC, [-]')
% ylabel('R_{0}, [Ohm]')
% title('Internal Resistance vs. SOC')
% plot(SOC_new, R0_new, 'LineWidth', 2)
% ylim([0 0.06])
% legend('Cell 1', 'Cell 3', 'Cell 6')
% box on
% set(gca, "FontSize", 24)


% Average to get the new R0(SOC) curve
R0_avg = mean(R0_new,2);

% Use 6th-order polynomial function:
B3 = polyfit(SOC_new, R0_avg, 6);
F3 = polyval(B3, SOC_new);
RMSE3 = sqrt(mean((R0_avg' - F3).^2));
disp(RMSE3)
% Results:
% Linear model Poly6:
% val(x) = p1*x^6 + p2*x^5 + p3*x^4 + p4*x^3 + p5*x^2 + p6*x + p7
% Coefficients (with 95% confidence bounds):
% p1 =       1.232  (1.124, 1.34)
% p2 =      -4.338  (-4.693, -3.983)
% p3 =       6.164  (5.706, 6.622)
% p4 =      -4.516  (-4.808, -4.223)
% p5 =       1.805  (1.709, 1.901)
% p6 =     -0.3779  (-0.3929, -0.3629)
% p7 =     0.06567  (0.06482, 0.06652)
% RMSE: 4.8034e-05
% for SOC in [0.0949, 1]


% Plot
hFig = figure(5);
set(hFig, 'Position', [100 100 600 500])
box on
hold on
% plot(SOC_new, R0_avg, '.k', 'MarkerSize', 10)
plot(SOC_new, R0_new, '.', 'MarkerSize', 15)
plot(SOC_new, F3, 'k', 'LineWidth', 2.5)
xlabel('SOC [-]')
ylabel('$R_{0}$ [$\Omega$]')
ylim([0 0.06])
legend({'Cell 1', 'Cell 3', 'Cell 6', 'Averaged Polynomial Fit'}, 'Location', 'best', 'interpreter', 'latex')
% title('Internal Resistance vs. SOC')
box on
set(gca, "FontSize", 28)

%% Test voc, soc functions
SOC_test = linspace(0,1);
for i = 1:length(SOC_test)
    [v_oc(i),dv_oc(i)] = voc(SOC_test(i));
end

for i = 1:length(SOC_test)
    [r_0(i),dr_0(i)] = r0(SOC_test(i));
end


hFig = figure(6);
set(hFig, 'Position', [100 100 900 500])
hold on
yyaxis left
plot(SOC_new, R0_new, '.', 'MarkerSize', 15)
plot(SOC_test, r_0, '-.', 'LineWidth', 2)
xlabel('SOC [-]')
ylabel('Internal Resistance [Ohm]')
yyaxis right
plot(SOC_test, dr_0, '-.', 'LineWidth', 2)
ylabel('dR0/dSOC [Ohm]')
legend('Cell 1', 'Cell 3', 'Cell 6', 'Avg Poly Fit', 'dR0/dSOC','Location', 'bestoutside')
set(gca, "FontSize", 28)

hFig = figure(7);
set(hFig, 'Position', [100 100 900 500])
hold on
yyaxis left
plot(SOC_new, OCV_new, '.', 'MarkerSize', 15)
plot(SOC_test, v_oc, '-.', 'LineWidth', 2)
xlabel('SOC [-]')
ylabel('Open Circuit Voltage [V]')
yyaxis right
plot(SOC_test, dv_oc, '-.', 'LineWidth', 2)
ylabel('dV_{OC}/dSOC [V]')
legend('Cell 1', 'Cell 3', 'Cell 6', 'Avg Poly Fit', 'dVoc/dSOC','Location', 'bestoutside')
set(gca, "FontSize", 28)











