%% Configuration file to load in all microgrid data
% Kevin Moy, 9/16/2021

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
% Control is DC power, bounded by kWAC rating of the LIB inverter

% Cell nominal voltage:
v_cell_nom = 3.63; % Volts

% Cell capacity:
Q_nom = 4.85; % Ampere-hours

% Number of cells in series
s = floor(LIB_INV_MAX_DC_V/v_cell_nom); %=228

% Number of cells in parallel
p = floor((LIB_NOM_ENG*1000)/(s*v_cell_nom*Q_nom)); %=118

% State limits
x_min = 0.2;
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
% yr_dt_sec = (YEAR_START:seconds(1):datetime(2021,12,31,23,59,59))';

%% Obtain load and PV data:
load('pv_gen.mat');
load('load_cons.mat');

