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
% temperature_data = readtable("./../Data/TJ_Temperature.xlsx");
boron_data = readtable("./../Data/TJ_d11B_pH.xlsx");

% Segregate into prepeturbation and during/after perturbation
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);
background_d11B4_uncertainty = nanstd(background_data.d11B);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4 = nanmin(perturbation_data.d11B);

d11B4_change = background_d11B4-perturbation_d11B4;

%% Set up distributions and sample them
% Choose the number of statistical samples
number_of_samples = 10000;

% Set up a distribution for each variable
% Uncertainties should be tuned to represent 2sigma uncertainty
% If not known, be generous
initial_co2_distribution = Geochemistry_Helpers.Distribution(0:10e-6:1000e-6,"Gaussian",[500e-6,50e-6]).normalise();
initial_saturation_distribution = Geochemistry_Helpers.Distribution(6:0.1:15,"Gaussian",[10.7,0.15]).normalise();
initial_temperature_distribution = Geochemistry_Helpers.Distribution(0:0.01:30,"Gaussian",[14.8+3.3*2,1.65]).normalise();
temperature_change_distribution = Geochemistry_Helpers.Distribution(0:0.1:10,"Gaussian",[4.2-1.1*2,0.55]).normalise();
ca_distribution = Geochemistry_Helpers.Distribution(0:0.1:20,"Flat",[8,10]).normalise();
mg_distribution = Geochemistry_Helpers.Distribution(20:0.1:61,"Flat",[50,52]).normalise();
epsilon_distribution = Geochemistry_Helpers.Distribution(27:0.01:31,"Gaussian",[27.2+0.6*2,0.3]).normalise();
initial_d11B_distribution = Geochemistry_Helpers.Distribution(0:0.1:40,"Gaussian",[background_d11B4,background_d11B4_uncertainty]).normalise();

% Create samplers for those distributions
initial_co2_sampler = Geochemistry_Helpers.Sampler(initial_co2_distribution,"latin_hypercube");
initial_saturation_sampler = Geochemistry_Helpers.Sampler(initial_saturation_distribution,"latin_hypercube");
initial_temperature_sampler = Geochemistry_Helpers.Sampler(initial_temperature_distribution,"latin_hypercube");
temperature_change_sampler = Geochemistry_Helpers.Sampler(temperature_change_distribution,"latin_hypercube");
ca_sampler = Geochemistry_Helpers.Sampler(ca_distribution,"latin_hypercube");
mg_sampler = Geochemistry_Helpers.Sampler(mg_distribution,"latin_hypercube");
epsilon_sampler = Geochemistry_Helpers.Sampler(epsilon_distribution,"latin_hypercube");
initial_d11B_sampler = Geochemistry_Helpers.Sampler(initial_d11B_distribution,"latin_hypercube");

% Get samples
initial_co2_sampler.getSamples(number_of_samples).shuffle();
initial_saturation_sampler.getSamples(number_of_samples).shuffle();
initial_temperature_sampler.getSamples(number_of_samples).shuffle();
temperature_change_sampler.getSamples(number_of_samples).shuffle();
ca_sampler.getSamples(number_of_samples).shuffle();
mg_sampler.getSamples(number_of_samples).shuffle();
epsilon_sampler.getSamples(number_of_samples).shuffle();
initial_d11B_sampler.getSamples(number_of_samples).shuffle();

% Make sure we don't pick 0 CO2 - could also find a better distribution
% than Gaussian so this doesn't happen
initial_co2_sampler.samples(initial_co2_sampler.samples==0)=60e-6;

%% Process background
% Set up an array to do the d11B->CO2 calculations
initial = BuCC.d11BCO2().create(number_of_samples);

% Assign all the necessary values
initial.collate("species_calibration").collate("d11B_measured").assignToEach("value",initial_d11B_sampler.samples);

initial.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial.collate("boron").collate("d11B_sw").assignToAll("value",NaN);
initial.collate("boron").assignToEach("epsilon",epsilon_sampler.samples);

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

%% d11B_sw
initial_d11B_sw = initial.collate("boron").collate("d11B_sw").collate("value");
initial_d11B_sw_distribution = Geochemistry_Helpers.Distribution.fromSamples((0:0.1:40)',initial_d11B_sw).normalise();

% We also have a d11B_sw constraint from the range of measured d11B
% So get the minimum and maximum
raw_d11B_minimum = min(boron_data.d11B);
raw_d11B_maximum = max(boron_data.d11B);

% And their uncertainties
raw_d11B_minimum_uncertainty = boron_data.uncertainty(boron_data.d11B==raw_d11B_minimum);
raw_d11B_maximum_uncertainty = boron_data.uncertainty(boron_data.d11B==raw_d11B_maximum);

% Create distributions and samplers
raw_d11B_minimum_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"Gaussian",[raw_d11B_minimum,raw_d11B_minimum_uncertainty]).normalise();
raw_d11B_maximum_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"Gaussian",[raw_d11B_maximum,raw_d11B_maximum_uncertainty]).normalise();
raw_d11B_minimum_sampler = Geochemistry_Helpers.Sampler(raw_d11B_minimum_distribution,"latin_hypercube");
raw_d11B_maximum_sampler = Geochemistry_Helpers.Sampler(raw_d11B_maximum_distribution,"latin_hypercube");

% Get the samples
raw_d11B_minimum_sampler.getSamples(number_of_samples).shuffle();
raw_d11B_maximum_sampler.getSamples(number_of_samples).shuffle();

% The d11B_sw can be anywhere between measured d11B and measured d11B +
% epsilon
% So get random values between 0 and 1 to represent what fraction of
% epsilon is expressed
fraction_of_epsilon_distribution = Geochemistry_Helpers.Distribution(-0.1:0.1:1.1,"flat",[0,1]).normalise();
fraction_of_epsilon_sampler = Geochemistry_Helpers.Sampler(fraction_of_epsilon_distribution,"latin_hypercube");
fraction_of_epsilon_sampler.getSamples(number_of_samples).shuffle();

% Calculate estimates of d11B_sw from the minimum
d11B_sw_from_minimum = raw_d11B_minimum_sampler.samples+(fraction_of_epsilon_sampler.samples.*epsilon_sampler.samples);

% Shuffle the epsilons so we don't get paired samples
fraction_of_epsilon_sampler.shuffle();
% Calculate estimates of d11B_sw from the maximum
d11B_sw_from_maximum = raw_d11B_maximum_sampler.samples+(fraction_of_epsilon_sampler.samples.*epsilon_sampler.samples);

% Turn the samples into a distribution
d11B_sw_from_minimum_distribution = Geochemistry_Helpers.Distribution.fromSamples(initial_d11B_sw_distribution.bin_edges,d11B_sw_from_minimum).normalise();
d11B_sw_from_maximum_distribution = Geochemistry_Helpers.Distribution.fromSamples(initial_d11B_sw_distribution.bin_edges,d11B_sw_from_maximum).normalise();

% Multiply them together to get the mutually inclusive region - and put
% that into a distribution
d11B_sw_from_minimum_maximum = d11B_sw_from_minimum_distribution.probabilities.*d11B_sw_from_maximum_distribution.probabilities;
d11B_sw_from_minimum_maximum_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"manual",d11B_sw_from_minimum_maximum).normalise();

combined_initial_d11B_sw = d11B_sw_from_minimum_maximum_distribution.probabilities.*initial_d11B_sw_distribution.probabilities;
combined_initial_d11B_sw_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"manual",combined_initial_d11B_sw).normalise();


figure(1);
clf
hold on
d11B_sw_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:60,initial.collate("boron").collate("d11B_sw").collate("value")).normalise();

combined_initial_d11B_sw_distribution.plot();
d11B_sw_distribution.plot();
d11B_sw_from_minimum_maximum_distribution.plot("Color",[0.8,0.8,0]);

%% Process minimum d11B measured
% Set up an array to do calculations
after = BuCC.d11BCO2().create(number_of_samples);

% Assign all the necessary variables
after.collate("species_calibration").collate("d11B_measured").assignToEach("value",initial_d11B_sampler.samples-d11B4_change);

after.collate("boron").collate("pH").assignToAll("pValue",NaN);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));
after.collate("boron").assignToEach("epsilon",epsilon_sampler.samples);

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

pH_initial_99 = pH_initial_distribution.quantile(0.99);
pH_initial_95 = pH_initial_distribution.quantile(0.95);
pH_initial_50 = pH_initial_distribution.quantile(0.5);

d11Bsw_initial_50 = initial_d11B_sw_distribution.mean();
d11Bsw_initial_uncertainty = initial_d11B_sw_distribution.standard_deviation();

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
