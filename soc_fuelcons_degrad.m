%% Produce SOC vs. fuel consumption, SOC vs. degradation
% Kevin Moy, 9/16/2021
clearvars
close all
clc
set(0,'defaultTextInterpreter','latex');
%% Load in microgrid configuration
mg_config

%% Define initial, final times

t_ini_strs = {'09-Oct-2019 00:00:00', '09-Oct-2019 15:15:00', ...
            '09-Oct-2019 22:30:00', '10-Oct-2019 09:45:00'};

t_fin_strs = {'12-Oct-2019 17:45:00', '12-Oct-2019 10:15:00', ...
            '12-Oct-2019 12:30:00', '12-Oct-2019 05:30:00'};
        
% Define temperature for degradation model

temps = [13.56, 13.56, 13.56, 13.56];

%% Create load, PV vectors

for i = 1:4
    t_vec = [YEAR_START; datetime([t_ini_strs{i}; t_fin_strs{i}])];

    t_ind = ((minutes(t_vec(2:3) - t_vec(1)))/15) + 1;
    
    loads{i} = ld(t_ind(1):t_ind(2));
    pvs{i} = pv(t_ind(1):t_ind(2));
    dts{i} = yr_dt(t_ind(1):t_ind(2));
    dts_soc{i} = yr_dt(t_ind(1):t_ind(2)+1);
    outage_lens(i) = length(ld(t_ind(1):t_ind(2)));
end

%% SOCs to try
SOCs = linspace(x_min, x_max, 27);

% %% Store vector of initial costate values
% lambda_inits = zeros(length(SOCs),1);
% SOC_fins = zeros(length(SOCs),1);

%% Produce plots of degradation vs. fuel consumption for each outage
lambda_min = -150;
lambda_max = -80;
for k = 1:4
    % Retrieve outage load, PV generation, length
    ld1 = loads{k};
    pv1 = pvs{k};
    outage_len = outage_lens(k);
    
    % Reset all loop storage vectors
    SOC_fins = zeros(length(SOCs), 1);
    lambda_inits = zeros(length(SOCs), 1);
    u_opts = zeros(outage_len, length(SOCs));
    SOC_opts = zeros(outage_len+1, length(SOCs));
    fuel_consumps = zeros(length(SOCs), 1);
    cap_losses = zeros(length(SOCs), 1);
    
    % Search for optimal initial costate values
    disp('Finding optimal initial costate values for each value of SOC:')
    for i = 1:length(SOCs)
        disp(['Iteration ', num2str(i)])
        [s_init_vec, SOC_end_vec] = init_costate_search(outage_len, pv1, ld1, x_max, x_min, ...
                    LIB_INV_SIZE_KW, LIB_EFF_CHG, ...
                    socdot, P_batt_range_socdot, SOC_range_socdot, ...
                    dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
                    h, lambda_min, lambda_max, SOCs(i));

        SOC_fins(i) = SOC_end_vec(end);
        lambda_inits(i) = s_init_vec(end);
        disp(['SOC ', num2str(SOCs(i)), ' Optimal costate value ', num2str(lambda_inits(i)), ' found'])
    end
    disp('Done finding optimal initial costate values!')

    disp('Finding optimal fuel consumption, SOC, LIB power from optimal intial costate values:')
    for i = 1:length(lambda_inits)
        disp(['Iteration ', num2str(i)])
        [u_opt, x, lambda] = mgpmpecm(outage_len, pv1, ld1, x_max, x_min, LIB_INV_SIZE_KW, LIB_EFF_CHG, ...
        socdot, P_batt_range_socdot, SOC_range_socdot, ...
        dfdx, P_batt_range_dfdx, SOC_range_dfdx, ...
        h, lambda_inits(i));
        dg1 = max(ld1-pv1-u_opt,0);
        pv_curt = max(pv1-ld1-u_opt,0);
        x_f = x(end);

        u_opts(:,i) = u_opt;
        SOC_opts(:,i) = x;
        fuel_consumps(i) = sum(genset_model(dg1))*h;
    end
    disp('Done with that.')

    % Compute baseline fuel consumption

    L_base = sum(genset_model(max(ld1-pv1,0)))*h;
    disp('Finding degradation from optimal SOC trajectories:')
    for i = 1:length(SOCs)
        disp(['Iteration ', num2str(i)])
        [cap_loss] = emp_deg_model(SOC_opts(:,i), Q_nom, h, temps(k));
        cap_losses(:,i) = cap_loss;
    end
    disp('Done with that.')
    u_opts_ks{k} = u_opts;
    SOC_opts_ks{k} = SOC_opts;
    L_pcts{k} = fuel_consumps/L_base;
    cap_pcts{k} = cap_losses/Q_nom;
end

%% Plot Pareto Front
%indices to plot given restrictions on Andrea's model 
%(SOC_min = [0.2,0.45])

inds = find(abs(SOCs-0.2) < 0.001):find(abs(SOCs-0.45) < 0.001);

len = length(inds);
% red = [0.9, 0.9, 0.9];
% blue = [0.2, 0.2, 0.2];
red = [1, 0, 0];
blue = [0, 0, 1];
colors_p = [linspace(red(1),blue(1),len)', linspace(red(2),blue(2),len)', linspace(red(3),blue(3),len)'];

markers = {'s', 'o', 'd', '^'};

% hFig = figure(1);
% set(hFig, 'Position', [100 100 600 500])
% hold on
% for j = 1:4
%     fuel_pct = L_pcts{j};
%     cap_pct = cap_pcts{j};
% %     plot(100-fuel_pct(inds)*100, cap_pct(1,inds)*100, markers{j}, 'MarkerSize', 10)
% %     plot(100-fuel_pct(inds)*100, cap_pct(1,inds)*100, 'LineWidth', 2)
%     scatter(100-fuel_pct(inds)*100, cap_pct(1,inds)*100, 140, markers{j}, 'filled')
% end
% % title({"Total Fuel Consumption and", "Cell Capacity Loss vs. Final LIB SOC"})
% xlabel("Reduction in Diesel consumption [\%]")
% ylabel("Reduction in Cell Capacity [\%]")
% ylim([1, 2.2])
% % plot(lambda_opts, 'LineWidth', 2)
% [l, hobj, hout, mout] = legend('Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'location', 'east', 'interpreter', 'latex', 'FontSize', 28);
% M = findobj(hobj,'type','patch');
% set(M,'MarkerSize',sqrt(140));
% l.FontSize = 28;
% set(gca, "FontSize", 28)
% set(gca,'TickLabelInterpreter','latex')
% box on


% hFig = figure(3);
% set(hFig, 'Position', [100 100 600 500])
% hold on
% for j = 1:4
%     fuel_pct = L_pcts{j};
%     cap_pct = cap_pcts{j};
%     plot(100-fuel_pct(inds)*100, cap_pct(1,inds)*100, 'LineWidth', 2)
% end
% % title({"Total Fuel Consumption and", "Cell Capacity Loss vs. Final LIB SOC"})
% xlabel("Reduction in Diesel Consumption [\%]")
% ylabel("Reduction in Cell Capacity [\%]")
% ylim([1, 2.2])
% % plot(lambda_opts, 'LineWidth', 2)
% legend('Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'location', 'east', 'interpreter', 'latex')
% set(gca, "FontSize", 28)

% hFig = figure(4);
% set(hFig, 'Position', [100 100 600 500])
% hold on
% for j = 1
%     fuel_pct = L_pcts{j};
%     cap_pct = cap_pcts{j};
%     for i = inds
%         plot(100-fuel_pct(i)*100, cap_pct(1,i)*100, markers{j}, 'Color', colors_p(i,:), 'MarkerSize', 10, 'MarkerFaceColor', colors_p(i,:), 'LineWidth', 1.5)
%     end
% %     plot(100-fuel_pct(inds)*100, cap_pct(1,inds)*100, 'k', 'LineWidth', 2)
% end
% % title({"Total Fuel Consumption", "and Cell Capacity Loss", "vs. Final LIB SOC, Phase 1"})
% ylim([1, 2.2])
% xlabel("Reduction in Diesel Consumption [\%]")
% ylabel("Reduction in Cell Capacity [\%]")
% % plot(lambda_opts, 'LineWidth', 2)
% legend('$SOC(t_f) = 20\%$', '','','','','','','','','','$SOC(t_f)=45\%$', 'location', 'east', 'interpreter', 'latex')
% % legend('$SOC(t_f) = 20\%$', '','','','','','','','','','','','', ...
% %     '','','','','','','','','','','','','','$SOC(t_f) = 85\%$', 'location', 'best', 'interpreter', 'latex')
% set(gca, "FontSize", 28)
% set(gca,'TickLabelInterpreter','latex')
% box on

hFig = figure(5);
set(hFig, 'Position', [100 100 600 500])
hold on
for j = 1:4
    fuel_pct = L_pcts{j};
    cap_pct = cap_pcts{j};
    for i = inds
        plot(100-fuel_pct(i)*100, cap_pct(1,i)*100, markers{j}, 'Color', colors_p(i,:), 'MarkerSize', 15, 'MarkerFaceColor', colors_p(i,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5)
    end
%     plot(100-fuel_pct(inds)*100, cap_pct(1,inds)*100, 'k', 'LineWidth', 2)
end
% title({"Total Fuel Consumption", "and Cell Capacity Loss", "vs. Final LIB SOC, Phase 1"})
ylim([1, 2.2])
xlabel("Reduction in Diesel Consumption [\%]")
ylabel("Reduction in Cell Capacity [\%]")
% plot(lambda_opts, 'LineWidth', 2)
% legend({"Phase 1" + newline + "$SOC(t_f) =0.2$", '','','','','','','','','',"$SOC(t_f) =0.45$", ...
%     "Phase 3" + newline + "$SOC(t_f) =0.2$", '','','','','','','','','',"$SOC(t_f) =0.45$"}, 'location', 'east', 'interpreter', 'latex')
legend({"Phase 1", '','','','','','','','','','', ...
    "Phase 2", '','','','','','','','','','', ...
    "Phase 3", '','','','','','','','','','', ...
    "Phase 4", '','','','','','','','','',""}, 'location', 'eastoutside', 'interpreter', 'latex')
% legend('$SOC(t_f) = 20\%$', '','','','','','','','','','','','', ...
%     '','','','','','','','','','','','','','$SOC(t_f) = 85\%$', 'location', 'best', 'interpreter', 'latex')
set(gca, "FontSize", 26)
set(gca,'TickLabelInterpreter','latex')
box on



%% Plot for particular Phase 1 SOC evolutions

len = length(inds);
% red = [0.7, 0.7, 0.7];
% blue = [0.1, 0.1, 0.1];
red = [1, 0, 0];
blue = [0, 0, 1];
colors_p = [linspace(red(1),blue(1),len)', linspace(red(2),blue(2),len)', linspace(red(3),blue(3),len)'];


soc_ph_1 = SOC_opts_ks{1};
hFig = figure(6);
set(hFig, 'Position', [100 100 600 500])
hold on
% socs_plot_ind = [1 9 13 17 21 27];
socs_plot_ind = [1 5 7 11];
for i = socs_plot_ind
    plot(dts_soc{1},soc_ph_1(:,i), 'Color', colors_p(i,:), 'LineWidth', 2.5)
end
% plot(dts_soc{1}, x_max*ones(length(dts_soc{1})), '--k', 'LineWidth', 2)
% plot(dts_soc{1}, x_min*ones(length(dts_soc{1})), ':k', 'LineWidth', 2)
% title({"SOC Trajectory", "vs. Final LIB SOC, Phase 1"})
ylabel("SOC [-]")
% ylim([0 1])
ylim([0.2 0.85])
% plot(lambda_opts, 'LineWidth', 2)
legendStrings = [string(SOCs(socs_plot_ind))];
% legendStrings = [string(SOCs(socs_plot_ind)), '$SOC_{UB}$', '$SOC_{LB}$'];
lgd = legend(legendStrings, 'location', 'eastoutside', 'interpreter', 'latex');
lgd.Title.String = '$SOC_f$';
set(gca, "FontSize", 28)
set(gca,'TickLabelInterpreter','latex')
box on



