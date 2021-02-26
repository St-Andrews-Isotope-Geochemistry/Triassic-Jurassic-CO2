% Minimum pH change sampler
% There are uncertainties on a lot of parameters necessary to calculate the
% change in pH
% To account for those uncertainties we want distributions for:
% - initial temperature
% - initial saturation state
% - initial CO2
% - ca
% - mg
% - temperature change
% And to set a 
% - worst case scenario epsilon (i.e. for minimum pH change)

% We'll work at 95% uncertainty (2sigma)
clear

%% Load in data
temperature_data = readtable("./../Data/TJ_Temperature.xlsx");
boron_data = readtable("./../Data/TJ_d11B_pH.xlsx");

% Segregate into prepeturbation and during/after perturbation
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4 = nanmin(perturbation_data.d11B);

%% Set up distributions and sample them
% Choose the number of statistical samples
number_of_samples = 10000;

% Set up a distribution for each variable
% Uncertainties should be tuned to represent 2sigma uncertainty
% If not known, be generous
initial_temperature_distribution = Geochemistry_Helpers.Distribution(0:2:40,"Gaussian",[20,10]).normalise();
initial_saturation_distribution = Geochemistry_Helpers.Distribution(6:1:15,"Gaussian",[10.7,1]).normalise();
initial_co2_distribution = Geochemistry_Helpers.Distribution(0:100e-6:1000e-6,"Gaussian",[500e-6,200e-6]).normalise();
ca_distribution = Geochemistry_Helpers.Distribution(0:0.5:20,"Flat",[1,17]).normalise();
mg_distribution = Geochemistry_Helpers.Distribution(26:0.5:50,"Flat",[28,42]).normalise();
temperature_change_distribution = Geochemistry_Helpers.Distribution(-20:5:30,"Gaussian",[10,5]).normalise();

% Create samplers for those distributions
initial_temperature_sampler = Geochemistry_Helpers.Sampler(initial_temperature_distribution,"latin_hypercube");
initial_saturation_sampler = Geochemistry_Helpers.Sampler(initial_saturation_distribution,"latin_hypercube");
initial_co2_sampler = Geochemistry_Helpers.Sampler(initial_co2_distribution,"latin_hypercube");
ca_sampler = Geochemistry_Helpers.Sampler(ca_distribution,"latin_hypercube");
mg_sampler = Geochemistry_Helpers.Sampler(mg_distribution,"latin_hypercube");
temperature_change_sampler = Geochemistry_Helpers.Sampler(temperature_change_distribution,"latin_hypercube");

% Get samples
initial_temperature_sampler.getSamples(number_of_samples).shuffle();
initial_saturation_sampler.getSamples(number_of_samples).shuffle();
initial_co2_sampler.getSamples(number_of_samples).shuffle();
ca_sampler.getSamples(number_of_samples).shuffle();
mg_sampler.getSamples(number_of_samples).shuffle();
temperature_change_sampler.getSamples(number_of_samples).shuffle();

% Make sure we don't pick 0 CO2 - could also find a better distribution
% that Gaussian so this doesn't happen
initial_co2_sampler.samples(initial_co2_sampler.samples==0)=60e-6;

%% Process background
% Set up an array to do the d11B->CO2 calculations
initial = BuCC.d11BCO2().create(number_of_samples);

% Assign all the necessary values
initial.collate("species_calibration").collate("d11B_measured").assignToAll("value",background_d11B4);

initial.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial.collate("boron").assignToAll("epsilon",27.2-(0.6*2));
initial.collate("boron").collate("d11B_sw").assignToAll("value",NaN);

cc = initial.collate("carbonate_chemistry");
cc.assignToAll("units"," mol/kg");
cc.assignToEach("temperature",initial_temperature_sampler.samples);
cc.assignToAll("salinity",35);
cc.assignToAll("oceanic_pressure",0);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToEach("calcium",ca_sampler.samples);
cc.assignToEach("magnesium",mg_sampler.samples);

initial.collate("carbonate_chemistry").collate("atmospheric_co2").assignToEach("partial_pressure",initial_co2_sampler.samples);
initial.collate("carbonate_chemistry").assignToEach("saturation_state",initial_saturation_sampler.samples);

myami = MyAMI.MyAMI("Precalculated",true);
initial.collate("carbonate_chemistry").collate("equilibrium_coefficients").assignToAll("MyAMI",myami);

% Do carbonate chemistry calculations
initial.calculate();

%% Process minimum d11B measured
% Set up an array to do calculations
after = BuCC.d11BCO2().create(number_of_samples);

% Assign all the necessary variables
after.collate("species_calibration").collate("d11B_measured").assignToAll("value",perturbation_d11B4);

after.collate("boron").collate("pH").assignToAll("pValue",NaN);
after.collate("boron").assignToAll("epsilon",27.2-(0.6*2));
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("units"," mol/kg");
cc.assignToEach("temperature",initial_temperature_sampler.samples+temperature_change_sampler.samples);
cc.assignToAll("salinity",35);
cc.assignToAll("oceanic_pressure",0);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToEach("calcium",ca_sampler.samples);
cc.assignToEach("magnesium",mg_sampler.samples);

after.collate("carbonate_chemistry").assignToEach("alkalinity",initial.collate("carbonate_chemistry").collate("alkalinity"));
after.collate("carbonate_chemistry").collate("equilibrium_coefficients").assignToAll("MyAMI",myami);

% Do carbonate chemistry calculations
after.calculate();

%% Postprocessing
% Collate the information
pH_after = after.collate("carbonate_chemistry").collate("pH").collate("pValue");
pH_initial = initial.collate("carbonate_chemistry").collate("pH").collate("pValue");

% If a value in either array is imaginary then we don't want it
real_pH = imag(pH_after)==0 & imag(pH_initial)==0;

% Create distributions from the samples
pH_after_distribution = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,pH_after(real_pH)).normalise();
pH_initial_distribution = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,pH_initial(real_pH)).normalise();

pH_change = pH_after(real_pH) - pH_initial(real_pH);
pH_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(-1:0.01:1,pH_change).normalise();

%% Statistics
% Can also calculate metrics to describe the distribution
pH_99 = pH_change_distribution.quantile(0.99);
pH_95 = pH_change_distribution.quantile(0.95);
pH_50 = pH_change_distribution.quantile(0.5);

%% Plotting
figure(1)
clf
hold on
pH_initial_distribution.plot();
pH_after_distribution.plot();
legend(["Background","Perturbation"],'Location','NorthWest','Box','Off');

xlabel("pH");
ylabel("Likelihood");


figure(2)
clf
hold on
pH_change_distribution.plot();

xlabel("\Delta pH");
ylabel("Likelihood");
