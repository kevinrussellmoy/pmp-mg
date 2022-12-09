%% Convert SAM PV 1 hour data to 15-minute PV data
% Kevin Moy
% Sep 8 2021

clearvars
close all
clc

%% Load SAM data based on BLRM specs (p. 84 in CEC report)

pv_8760 = readtable('pv_8760_data.csv');

pv = repelem(table2array(pv_8760(:,2)), 4); % Replicate 4x for 15-minute data

% plot(pv_15_min)

save('pv_gen.mat', 'pv')
