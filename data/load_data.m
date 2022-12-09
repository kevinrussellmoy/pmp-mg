%% Load in load data and scale appropriately
% Kevin Moy
% Sep 8 2021
%% Load residential and commercial data\

resi = readtable('resi_load.csv');
comm = readtable('bldg_load.csv');

%% Sum up resi loads and multiply by 10 for 5*18 = 90 homes
resi_sum = sum(table2array(resi(:, 2:end)), 2)*5;
plot(resi_sum) 

%% Include scaled down hospital (0.25x)

hosp = table2array(comm(:,2))*0.25;
plot(hosp)

%% Sum together resi and load
ld = resi_sum + hosp;
disp(max(ld))
plot(ld)
% With resi of 5*18 = 90 homes, and 0.25 hospital, peak load is 350.2675 kW

save('load_cons.mat', 'ld')