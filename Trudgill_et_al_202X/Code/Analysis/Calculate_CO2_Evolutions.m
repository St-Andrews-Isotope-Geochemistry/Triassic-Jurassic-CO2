% TJ_CO2 v3
% Using: Initial CO2 between 500ppm and 5000ppm
%        Omega between 5 and 11
% Gives: Viable initial pH distribution
%        (Viable initial alkalinity distribution)
%
% Using: Initial pH distribution
%        Measured d11B
%        Epsilon + Kb
% Gives: d11B_sw distribution
%
% Using: Possible d11B_measured evolutions + species calibration
% Gives: Possible d11B_4 evolutions
%
% Using: d11B_4 range + epsilon
% Gives: d11B_sw distribution
%
% Combine both d11B_sw distributions but avoid resampling to retain
% internal consistency
%
% Using: d11B_sw + already chosen parameters
% Gives: pH
%
% Using: Assumptions about viable saturation state + pH
% Gives: Plausible alkalinity
%
% Using: Assumptions on shape and rate of change of alkalinity
% Gives: Possible alkalinity evolutions
%        
% Combine both alkalinity possibilities (again avoid resampling)
%
% Using: Alkalinity evolutions + pH
% Gives: CO2 evolutions

clear

%%
tic
data_directory = "./../../Data/";

% Choose the number of statistical samples
number_of_samples = 1000;

% Load in the data
boron_data = readtable(data_directory+"TJ_d11B.xlsx","Sheet","Delta_Temperature");
raw_input_parameters = fileread(data_directory+"CO2_Evolutions_Input.json");
input_parameters = jsondecode(raw_input_parameters);

raw_alkalinity_constraints = fileread(data_directory+"Alkalinity_Constraints.json");
alkalinity_constraints = jsondecode(raw_alkalinity_constraints);

for alkalinity_constraint_index = 1:size(alkalinity_constraints,1)
    alkalinity_constraint_distributions(alkalinity_constraint_index) = Geochemistry_Helpers.Distribution(-1:0.0001:1,"Gaussian",alkalinity_constraints(alkalinity_constraint_index,2:3)).normalise();
    alkalinity_constraint_distributions(alkalinity_constraint_index).location = alkalinity_constraints(alkalinity_constraint_index,1);
end


%% Distributions and Samplers
% Generate distributions for each variable to represent uncertainties
initial_co2_distribution = Geochemistry_Helpers.Distribution(0:100e-6:10000e-6,"Flat",input_parameters.co2/1e6).normalise();
initial_omega_distribution = Geochemistry_Helpers.Distribution(0:0.1:12,"Flat",input_parameters.saturation_state).normalise();
ca_distribution = Geochemistry_Helpers.Distribution(0:0.1:20,"Flat",input_parameters.calcium).normalise();
mg_distribution = Geochemistry_Helpers.Distribution(20:0.1:61,"Flat",input_parameters.magnesium).normalise();
species_calibration_gradient_distribution = Geochemistry_Helpers.Distribution(0:0.01:1.1,"Flat",[0.1,1]).normalise();
species_calibration_offset_distribution = Geochemistry_Helpers.Distribution(-10:0.1:10,"Flat",[-5,5]).normalise();

initial_temperature_distributions = [Geochemistry_Helpers.Distribution(-10:0.1:30,"Gaussian",input_parameters.initial_temperature(1,2:3)).normalise(input_parameters.initial_temperature(1,1)),Geochemistry_Helpers.Distribution(-10:0.1:30,"Gaussian",input_parameters.initial_temperature(2,2:3)).normalise(input_parameters.initial_temperature(2,1))];
initial_temperature_distribution = Geochemistry_Helpers.Distribution(-10:0.1:30,"Manual",initial_temperature_distributions(1).probabilities+initial_temperature_distributions(2).probabilities).normalise();

epsilon_distributions = [Geochemistry_Helpers.Distribution(23:0.01:30,"Gaussian",input_parameters.epsilon(1,2:3)).normalise(input_parameters.epsilon(1,1)), Geochemistry_Helpers.Distribution(23:0.01:30,"Gaussian",input_parameters.epsilon(2,2:3)).normalise(input_parameters.epsilon(2,1))];
epsilon_distribution = Geochemistry_Helpers.Distribution(23:0.01:30,"Manual",epsilon_distributions(1).probabilities+epsilon_distributions(2).probabilities).normalise();

d11B_measured_distributions = Geochemistry_Helpers.Distribution().create([height(boron_data),1]);
for d11B_distribution_index = 1:numel(d11B_measured_distributions)
    d11B_measured_distributions(d11B_distribution_index) = Geochemistry_Helpers.Distribution(0:0.1:40,"Gaussian",[boron_data.d11B(d11B_distribution_index),boron_data.d11B_uncertainty(d11B_distribution_index)/2]).normalise();
end
temperature_distributions = Geochemistry_Helpers.Distribution().create([height(boron_data),1]);
for temperature_distribution_index = 1:numel(d11B_measured_distributions)
    temperature_distributions(temperature_distribution_index) = Geochemistry_Helpers.Distribution(-20:0.1:20,"Gaussian",[boron_data.delta_temperature(temperature_distribution_index),boron_data.delta_temperature_uncertainty(temperature_distribution_index)]).normalise();
end

% Use a sampling technique for those distributions
initial_co2_sampler = Geochemistry_Helpers.Sampler(initial_co2_distribution,"latin_hypercube");
initial_omega_sampler = Geochemistry_Helpers.Sampler(initial_omega_distribution,"latin_hypercube");
initial_temperature_sampler = Geochemistry_Helpers.Sampler(initial_temperature_distribution,"latin_hypercube");
ca_sampler = Geochemistry_Helpers.Sampler(ca_distribution,"latin_hypercube");
mg_sampler = Geochemistry_Helpers.Sampler(mg_distribution,"latin_hypercube");
epsilon_sampler = Geochemistry_Helpers.Sampler(epsilon_distribution,"latin_hypercube");
species_calibration_gradient_sampler = Geochemistry_Helpers.Sampler(species_calibration_gradient_distribution,"latin_hypercube");
species_calibration_offset_sampler = Geochemistry_Helpers.Sampler(species_calibration_offset_distribution,"latin_hypercube");

d11B_measured_samplers = Geochemistry_Helpers.Sampler().create([height(boron_data),1]);
for d11B_sampler_index = 1:numel(d11B_measured_distributions)
    d11B_measured_samplers(d11B_sampler_index) = Geochemistry_Helpers.Sampler(d11B_measured_distributions(d11B_sampler_index),"latin_hypercube");
end
temperature_samplers = Geochemistry_Helpers.Sampler().create([height(boron_data),1]);
for temperature_sampler_index = 1:numel(temperature_distributions)
    temperature_samplers(temperature_sampler_index) = Geochemistry_Helpers.Sampler(temperature_distributions(temperature_sampler_index),"latin_hypercube");
end

% Get the samples
initial_co2_sampler.getSamples(number_of_samples).shuffle();
initial_omega_sampler.getSamples(number_of_samples).shuffle();
initial_temperature_sampler.getSamples(number_of_samples).shuffle();
ca_sampler.getSamples(number_of_samples).shuffle();
mg_sampler.getSamples(number_of_samples).shuffle();
epsilon_sampler.getSamples(number_of_samples).shuffle();
species_calibration_gradient_sampler.getSamples(number_of_samples).shuffle();
species_calibration_offset_sampler.getSamples(number_of_samples).shuffle();

for d11B_sampler_index = 1:numel(d11B_measured_samplers)
    d11B_measured_samplers(d11B_sampler_index).getSamples(number_of_samples).shuffle();
end
d11B_measured_evolutions = d11B_measured_samplers.collate("samples");
for temperature_sampler_index = 1:numel(temperature_samplers)
    temperature_samplers(temperature_sampler_index).getSamples(number_of_samples).shuffle();
end
temperature_evolutions = temperature_samplers.collate("samples");
temperature_evolutions(1,:) = 0;

% Additional parameters
oceanic_pressure = input_parameters.oceanic_pressure;
atmospheric_pressure = input_parameters.atmospheric_pressure;
salinity = input_parameters.salinity;

%% Initial pH
initial_d11B_CO2 = BuCC.d11BCO2().create(number_of_samples);

% Species calibration
initial_d11B_CO2.species_calibration.assignToEach("coefficients",[species_calibration_gradient_sampler.samples;species_calibration_offset_sampler.samples]');
initial_d11B_CO2.species_calibration.d11B_measured.assignToEach("value",d11B_measured_samplers(1).samples);

% Boron
initial_d11B_CO2.boron.pH.assignToAll("pValue",NaN);
initial_d11B_CO2.boron.d11B_sw.assignToAll("value",NaN);
initial_d11B_CO2.boron.assignToEach("epsilon",epsilon_sampler.samples);

% Carbonate chemistry
% Main parameters
initial_d11B_CO2.carbonate_chemistry.atmospheric_co2.assignToEach("partial_pressure",initial_co2_sampler.samples);
initial_d11B_CO2.carbonate_chemistry.assignToEach("saturation_state",initial_omega_sampler.samples);
initial_d11B_CO2.carbonate_chemistry.assignToAll("units"," mol/kg");

% Ancillary
% Variable
initial_d11B_CO2.carbonate_chemistry.assignToEach("temperature",initial_temperature_sampler.samples);
initial_d11B_CO2.carbonate_chemistry.assignToEach("calcium",ca_sampler.samples);
initial_d11B_CO2.carbonate_chemistry.assignToEach("magnesium",mg_sampler.samples);

% Constant
initial_d11B_CO2.carbonate_chemistry.assignToAll("salinity",salinity);
initial_d11B_CO2.carbonate_chemistry.assignToAll("oceanic_pressure",oceanic_pressure);
initial_d11B_CO2.carbonate_chemistry.assignToAll("atmospheric_pressure",atmospheric_pressure);

% Create a MyAMI object
myami = MyAMI.MyAMI("Precalculated",true);
initial_d11B_CO2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

% Do carbonate chemistry calculations
initial_d11B_CO2.calculate();

% Collate the results of those calculations back to normal arrays
initial_pH = initial_d11B_CO2.carbonate_chemistry.pH.collate("pValue");
initial_alkalinity = initial_d11B_CO2.carbonate_chemistry.alkalinity;

%% We have starting conditions but can refine d11B_sw
d11B_4_evolutions = BuCC.BoronSpeciesCalibration().create([height(boron_data),number_of_samples]);
for evolution_index = 1:size(d11B_4_evolutions,2)    
    % Species calibration
    d11B_4_evolutions(:,evolution_index).assignToAll("form","polynomial");
    d11B_4_evolutions(:,evolution_index).assignToAll("coefficients",[species_calibration_gradient_sampler.samples(evolution_index);species_calibration_offset_sampler.samples(evolution_index)]');
    d11B_4_evolutions(:,evolution_index).d11B_measured.assignToEach("value",d11B_measured_evolutions(:,evolution_index));
    
    % Do species calculations
    d11B_4_evolutions(:,evolution_index).calculate();
end
d11B_4_evolutions_value = d11B_4_evolutions.d11B_4.value;

% Collate the viable d11B_sw's and make a distribution
initial_d11B_sw = initial_d11B_CO2.boron.d11B_sw.value;
initial_d11B_sw_sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.1:60,initial_d11B_sw,"latin_hypercube");

% We also have a d11B_sw constraint from the range of measured d11B
% So get the minimum and maximum
d11B_sw_minimum = max(d11B_4_evolutions_value);
d11B_sw_maximum = min(d11B_4_evolutions_value)+epsilon_sampler.samples;

%%
number_of_subsamples = number_of_samples;

d11B_sw_weighting = double(initial_d11B_sw>=d11B_sw_minimum' &  initial_d11B_sw<=d11B_sw_maximum');
[d11B_sw_resamples,resample_indices] = initial_d11B_sw_sampler.resample(d11B_sw_weighting,number_of_subsamples);

d11B_resamples = d11B_measured_evolutions(:,resample_indices);
initial_co2_resamples = initial_co2_sampler.samples(resample_indices);
initial_omega_resamples = initial_omega_sampler.samples(resample_indices);
initial_temperature_resamples = initial_temperature_sampler.samples(resample_indices);
ca_resamples = ca_sampler.samples(resample_indices);
mg_resamples = mg_sampler.samples(resample_indices);
epsilon_resamples = epsilon_sampler.samples(resample_indices);
species_calibration_gradient_resamples = species_calibration_gradient_sampler.samples(resample_indices);
species_calibration_offset_resamples = species_calibration_offset_sampler.samples(resample_indices);
temperature_resamples = temperature_evolutions(:,resample_indices);
alkalinity_resamples = initial_alkalinity(resample_indices);

initial_d11B_sw_redistribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:60,d11B_sw_resamples);
initial_co2_redistribution = Geochemistry_Helpers.Distribution.fromSamples(0:100e-6:10000e-6,initial_co2_resamples);
initial_omega_redistribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:12,initial_omega_resamples);
initial_temperature_redistribution = Geochemistry_Helpers.Distribution.fromSamples(-10:0.1:30,initial_temperature_resamples);
ca_redistribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,ca_resamples);
mg_redistribution = Geochemistry_Helpers.Distribution.fromSamples(20:0.1:61,mg_resamples);
epsilon_redistribution = Geochemistry_Helpers.Distribution.fromSamples(23:0.01:30,epsilon_resamples);
species_calibration_gradient_redistribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:1.1,species_calibration_gradient_resamples);
species_calibration_offset_redistribution = Geochemistry_Helpers.Distribution.fromSamples(-10:0.1:10,species_calibration_offset_resamples);
temperature_redistribution = Geochemistry_Helpers.Distribution.fromSamples(-20:0.1:20,temperature_resamples);

%%
pH_evolutions = BuCC.d11BCO2().create([height(boron_data),number_of_subsamples,2]);
for evolution_index = 1:size(pH_evolutions,2)
    % Low omega
    % Species calibration
    pH_evolutions(:,evolution_index,1).species_calibration.assignToAll("coefficients",[species_calibration_gradient_resamples(evolution_index);species_calibration_offset_resamples(evolution_index)]');
    pH_evolutions(:,evolution_index,1).species_calibration.d11B_measured.assignToEach("value",d11B_resamples(:,evolution_index));
    
    % Boron
    pH_evolutions(:,evolution_index,1).boron.pH.assignToAll("pValue",NaN);
    pH_evolutions(:,evolution_index,1).boron.d11B_sw.assignToAll("value",d11B_sw_resamples(evolution_index));
    pH_evolutions(:,evolution_index,1).boron.assignToAll("epsilon",epsilon_resamples(evolution_index));
    
    % Carbonate chemistry
    % Main parameters
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",NaN);
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.assignToEach("saturation_state",[5*ones(9,1);1*ones(13,1)]);
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("units"," mol/kg");
    
    % Ancillary
    % Variable
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.assignToEach("temperature",initial_temperature_resamples(evolution_index)+temperature_resamples(:,evolution_index));
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("calcium",ca_resamples(evolution_index));
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("magnesium",mg_resamples(evolution_index));
    
    % Constant
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("salinity",salinity);
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("oceanic_pressure",oceanic_pressure);
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("atmospheric_pressure",atmospheric_pressure);
    
    % Create a MyAMI object
    myami = MyAMI.MyAMI("Precalculated",true);
    pH_evolutions(:,evolution_index,1).carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);
    
    % Do carbonate chemistry calculations
    pH_evolutions(:,evolution_index,1).calculate();
    
    % High omega
    % Species calibration
    pH_evolutions(:,evolution_index,2).species_calibration.assignToAll("coefficients",[species_calibration_gradient_resamples(evolution_index);species_calibration_offset_resamples(evolution_index)]');
    pH_evolutions(:,evolution_index,2).species_calibration.d11B_measured.assignToEach("value",d11B_resamples(:,evolution_index));
    
    % Boron
    pH_evolutions(:,evolution_index,2).boron.pH.assignToAll("pValue",NaN);
    pH_evolutions(:,evolution_index,2).boron.d11B_sw.assignToAll("value",d11B_sw_resamples(evolution_index));
    pH_evolutions(:,evolution_index,2).boron.assignToAll("epsilon",epsilon_resamples(evolution_index));
    
    % Carbonate chemistry
    % Main parameters
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",NaN);
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("saturation_state",12);
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("units"," mol/kg");
    
    % Ancillary
    % Variable
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.assignToEach("temperature",initial_temperature_resamples(evolution_index)+temperature_resamples(:,evolution_index));
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("calcium",ca_resamples(evolution_index));
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("magnesium",mg_resamples(evolution_index));
    
    % Constant
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("salinity",salinity);
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("oceanic_pressure",oceanic_pressure);
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("atmospheric_pressure",atmospheric_pressure);
    
    % Create a MyAMI object
    myami = MyAMI.MyAMI("Precalculated",true);
    pH_evolutions(:,evolution_index,2).carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);
    
    % Do carbonate chemistry calculations
    pH_evolutions(:,evolution_index,2).calculate();  
end

% Collate the results of those calculations back to normal arrays
pH_evolutions_value = pH_evolutions.carbonate_chemistry.pH.pValue;
pH_evolutions_low_omega = pH_evolutions_value(:,:,1);

alkalinity_boundaries = pH_evolutions.carbonate_chemistry.alkalinity;

%% To calculate evolutions we need to know how alkalinity may have changed
% Do this using Gaussian process

gp = Geochemistry_Helpers.GaussianProcess("rbf",boron_data.age');
gp.observations = alkalinity_constraint_distributions;

gp.runKernel([0.5,0.1]);
gp.getSamples(number_of_samples);

scaled_alkalinity_samples = (((gp.samples.*8000)+alkalinity_resamples*1e6)/1e6)';

viable_alkalinity_samples = all(scaled_alkalinity_samples>=alkalinity_boundaries(:,:,1)) & all(scaled_alkalinity_samples<=alkalinity_boundaries(:,:,2));

alkalinity_evolutions_viable = scaled_alkalinity_samples(:,viable_alkalinity_samples);
pH_evolutions_viable = pH_evolutions_low_omega(:,viable_alkalinity_samples);
d11B_evolutions_viable = d11B_resamples(:,viable_alkalinity_samples);
species_calibration_gradient_viable = species_calibration_gradient_resamples(viable_alkalinity_samples);
species_calibration_offset_viable = species_calibration_offset_resamples(viable_alkalinity_samples);

d11B_sw_viable = d11B_sw_resamples(viable_alkalinity_samples);
epsilon_viable = epsilon_resamples(viable_alkalinity_samples);
initial_temperature_viable = initial_temperature_resamples(viable_alkalinity_samples);
ca_viable = ca_resamples(viable_alkalinity_samples);
mg_viable = mg_resamples(viable_alkalinity_samples);
temperature_viable = temperature_resamples(:,viable_alkalinity_samples);

%%
co2_evolutions = BuCC.d11BCO2().create([height(boron_data),size(alkalinity_evolutions_viable,2)]);
for evolution_index = 1:size(alkalinity_evolutions_viable,2)
    % Species calibration
    co2_evolutions(:,evolution_index,1).species_calibration.assignToAll("coefficients",[species_calibration_gradient_viable(evolution_index),species_calibration_offset_viable(evolution_index)]');
    co2_evolutions(:,evolution_index,1).species_calibration.d11B_measured.assignToEach("value",d11B_evolutions_viable(:,evolution_index));
    
    % Boron
    co2_evolutions(:,evolution_index,1).boron.pH.assignToAll("pValue",NaN);
    co2_evolutions(:,evolution_index,1).boron.d11B_sw.assignToAll("value",d11B_sw_viable(evolution_index));
    co2_evolutions(:,evolution_index,1).boron.assignToAll("epsilon",epsilon_viable(evolution_index));
    
    % Carbonate chemistry
    % Main parameters
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",NaN);
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.assignToEach("alkalinity",alkalinity_evolutions_viable(:,evolution_index));
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("units"," mol/kg");
    
    % Ancillary
    % Variable
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.assignToEach("temperature",initial_temperature_viable(evolution_index)+temperature_viable(:,evolution_index));
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("calcium",ca_viable(evolution_index));
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("magnesium",mg_viable(evolution_index));
    
    % Constant
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("salinity",salinity);
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("oceanic_pressure",oceanic_pressure);
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("atmospheric_pressure",atmospheric_pressure);
    
    % Create a MyAMI object
    myami = MyAMI.MyAMI("Precalculated",true);
    co2_evolutions(:,evolution_index,1).carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);
    
    % Do carbonate chemistry calculations
    co2_evolutions(:,evolution_index).calculate();
end
co2_evolutions_value = co2_evolutions.carbonate_chemistry.atmospheric_co2.x;
saturation_state_output = co2_evolutions.carbonate_chemistry.saturation_state;
alkalinity_output = co2_evolutions.carbonate_chemistry.alkalinity;
pH_output = co2_evolutions.carbonate_chemistry.pH.pValue;
