%% Test of SOC_dot derivative wrt. SOC map 
% based on cell current partial derivatives
% Kevin Moy
% 9/1/21

%% Problem parameters
clearvars
Q_nom = 4.85; % Ah
s = 228;
p = 118;

numpts_soc = 101;

numpts_P = 51;

min_P = -300;
max_P = 300;
min_SOC = 0.095;
max_SOC = 1;

% Range of P_batt:
P_batt_range = linspace(min_P,max_P,numpts_P);

% Range of SOC:
SOC_range = linspace(min_SOC,max_SOC,numpts_soc);

% [X, Y] =  meshgrid(P_batt_range, SOC_range)

for i = 1:numpts_soc
    soc = SOC_range(i);
    for j = 1:numpts_P
        P_batt = P_batt_range(j);
        di_dvoc = didv(soc, P_batt, s, p);
        di_dr0 = didr(soc, P_batt, s, p);
        [~, dvoc_dsoc] = voc(soc);
        [~, dr0_dsoc] = r0(soc);
%         disp(di_dvoc)
%         disp(di_dr0)
        
        dsoc_dot_dsoc(i,j) = -(1/Q_nom) * (di_dvoc*dvoc_dsoc + di_dr0*dr0_dsoc);
    end
end

% Plot

% hFig = figure(1000);
% plot(P_batt_range, dsoc_dot_dsoc)
% xlabel('Pack Power [kW]')
% ylabel('$\partial \dot{SOC} / \partial SOC$', 'Interpreter','latex')

% levs = linspace(min(min(dsoc_dot_dsoc)), max(max(dsoc_dot_dsoc)), 20);

% levs = [min(min(dsoc_dot_dsoc)), max(max(dsoc_dot_dsoc))];

% levs = [-flip(logspace(0, log10(-min(min(dsoc_dot_dsoc))), 10)), 0, logspace(0, log10(-min(min(dsoc_dot_dsoc))), 10)];

% maxx = max(max(dsoc_dot_dsoc));
% 
% minx = min(min(dsoc_dot_dsoc));
% 
% test = [-(-minx)^(1/3):0.4:0 0:0.4:maxx^(1/3)];
% 
% levs = test.^3;

levs = [-0.3, -0.2, -0.1, -0.04, 0, .04, 0.1, 0.2, 0.35];


hFig = figure(1);
set(hFig, 'Position', [100 100 600 500])
[c, h] = contour(P_batt_range./(s*p).*1000, SOC_range, dsoc_dot_dsoc, levs, 'LineWidth', 2);
h.LevelList=round(h.LevelList,3);
h.LevelList = sort([h.LevelList 0]);
clabel(c, h, 'FontSize', 22, 'labelspacing', 150)
% contour(P_batt_range, SOC_range, dsoc_dot_dsoc, levs, 'LineWidth', 2, 'ShowText','on');
a = colorbar;
a.Label.Interpreter = 'latex';
a.Label.String = '$\mathrm{\partial\dot{SOC}/\partial SOC}$ [s$^{-1}$]';
a.Label.FontSize = 28;
xlabel('Cell Power [WDC]', 'Interpreter','latex')
ylabel('SOC [-]', 'Interpreter','latex')
% title('$\partial \dot{SOC} / \partial SOC$', 'Interpreter','latex')
ylim([0.095 1])
ax = gca;
set(ax,'FontSize',28);
ax.Position = ax.Position - [0 0 .03 0];
box on

% hFig = figure(2);
% set(hFig, 'Position', [700 100 800 500])
% surf(P_batt_range, SOC_range, dsoc_dot_dsoc)
% colorbar
% xlabel('Pack Power [kW]')
% ylabel('SOC [-]')
% ylim([0.095 1])
% zlabel('$\partial \dot{SOC} / \partial SOC$', 'Interpreter','latex')
% box on
% ax = gca;
% set(ax,'FontSize',28);
% cb = colorbar('Location','eastoutside');
% ax.Position = ax.Position - [0 0 .1 .1];
% cb.Position = cb.Position + [.1 0 0 0];

%% Try plotting soc dot by itself??

clearvars soc_dot

for i = 1:numpts_soc
    for j = 1:numpts_P
        P_batt = P_batt_range(j);
        P_cell = P_batt/(s*p)*1000; % CONVERT TO WATTS
        soc = SOC_range(i);
        [v_oc,~] = voc(soc);
        [r_0, ~] = r0(soc);
        soc_dot(i,j) = -1/Q_nom * ( (v_oc - sqrt(v_oc^2-4*r_0*P_cell))/(2*r_0) ); 

%         soc_dot(i,j) = -1/Q_nom * (v_oc/(2*r_0) - sqrt((v_oc/(2*r_0))^2 - P_cell/r_0)); 
    end
end
hFig = figure(10);
set(hFig, 'Position', [100 100 600 500])
box on
hold on
% surf(P_batt_range, SOC_range, soc_dot)
[c, h] = contour(P_batt_range./(s*p).*1000, SOC_range, soc_dot, 8, 'LineWidth', 2);
h.LevelList=round(h.LevelList,3);
h.LevelList = sort([h.LevelList 0]);
clabel(c, h, 'FontSize', 22, 'labelspacing', 300)
a = colorbar;
a.Label.Interpreter = 'latex';
a.Label.String = '$\mathrm{\dot{SOC}}$ [s$^{-1}$]';
a.Label.FontSize = 28;
xlabel('Cell Power [WDC]', 'Interpreter','latex')
ylabel('SOC [-]', 'Interpreter','latex')
% zlabel('$\dot{SOC}$', 'Interpreter','latex')
ylim([0.095 1])
% ylabel('R_{0}, [Ohm]')
% ylim([0 0.06])
% legend('Cell 1', 'Cell 3', 'Cell 6', 'Averaged Polynomial Fit', 'Location', 'best')
% title('Internal Resistance vs. SOC')
ax = gca;
set(ax,'FontSize',28);
ax.Position = ax.Position - [0 0 .03 0];

%%
hFig = figure(11);
set(hFig, 'Position', [100 100 600 500])
box on
hold on
plot(P_batt_range, soc_dot)
xlabel('Pack Power [kWDC]', 'Interpreter','latex')
ylabel('$\dot{SOC}$ [1/s]', 'Interpreter','latex')
xlim([-300 300])
ylim([-1 1])
box on
set(gca, "FontSize", 28)

hFig = figure(12);
set(hFig, 'Position', [100 100 600 500])
box on
hold on
plot(P_batt_range, dsoc_dot_dsoc)
xlabel('Pack Power [kWDC]', 'Interpreter','latex')
ylabel('$\partial\dot{SOC}/\partial SOC$ [1/s]', 'Interpreter','latex')
xlim([-300 300])
ylim([-1 1])
box on
set(gca, "FontSize", 28)
