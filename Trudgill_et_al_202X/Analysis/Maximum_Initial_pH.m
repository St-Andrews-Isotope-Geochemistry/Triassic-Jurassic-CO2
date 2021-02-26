% Maximum Initial pH
% The most basic processing to calculate the maximum initial pH
% Use saturation state = 10.7, CO2 = 500ppm, and other parameters as
% specified
% This is done is both csys and my BuCC package to ensure they're the same
% A final check is done by using the output of BuCC to check that the
% saturation state comes out the same

%%
clear

%% Set up the variables
saturation_state_maximum = 10.7; % From Ridgwell modelling
co2_minimum = 500; % From Witkowski

temperature = 21.4; % From bulk temperatures
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = 17; % Horita
mg = 28; % Horita

%% Using csys
% Run random numbers just to get kspc
flag = 15;
myami = MyAMI.MyAMI("MyAMI",false);
[~,csys_results_kspc] = fncsysKMgCaV2(flag,temperature,salinity,pressure,NaN,NaN,NaN,NaN,2300,2000,NaN,mg,ca,myami);

% Calculate CO3 and CO2
co3 = ((saturation_state_maximum.*csys_results_kspc.Kspc)/(ca/1e3))*1e6;
ocean_co2 = co2_minimum*(csys_results_kspc.Kh);

% Calculate pH
flag = 3;
[~,csys_results] = fncsysKMgCaV2(flag,temperature,salinity,pressure,NaN,ocean_co2,NaN,co3,NaN,NaN,NaN,mg,ca,myami);

maximum_pH_csys = csys_results.pH;

%% Using BuCC
% Assign each of the parameters
carbonate_chemistry = BuCC.CarbonateChemistry();

carbonate_chemistry.temperature = temperature;
carbonate_chemistry.salinity = salinity;
carbonate_chemistry.oceanic_pressure = pressure;
carbonate_chemistry.atmospheric_pressure = 1;
carbonate_chemistry.calcium = ca;
carbonate_chemistry.magnesium = mg;

carbonate_chemistry.atmospheric_co2.partial_pressure = co2_minimum*1e-6;
carbonate_chemistry.saturation_state = saturation_state_maximum;

carbonate_chemistry.equilibrium_coefficients.MyAMI = myami;

% Calculate the carbonate system for pH
carbonate_chemistry.calculate();

maximum_pH_cc = carbonate_chemistry.pH.pValue;

%% Using BuCC package with pH + CO2 as input to check output saturation state
carbonate_chemistry_check = BuCC.CarbonateChemistry();

carbonate_chemistry_check.temperature = temperature;
carbonate_chemistry_check.salinity = salinity;
carbonate_chemistry_check.oceanic_pressure = pressure;
carbonate_chemistry_check.atmospheric_pressure = 1;
carbonate_chemistry_check.calcium = ca;
carbonate_chemistry_check.magnesium = mg;

carbonate_chemistry_check.atmospheric_co2.partial_pressure = co2_minimum*1e-6;
carbonate_chemistry_check.pH.pValue = maximum_pH_cc;

carbonate_chemistry_check.equilibrium_coefficients.MyAMI = myami;

carbonate_chemistry_check.calculate();

%% Difference
pH_difference = abs(maximum_pH_csys-maximum_pH_cc);
