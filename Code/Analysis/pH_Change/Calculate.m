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

% clear

%%
data_directory = "./../../../Data/";

% Choose the number of statistical samples
round_1.number_of_samples = 500;

% Load in the data
data.boron = readtable(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx","Sheet","Temperature_Calibrations");
data.boron = data.boron((~isnan(data.boron.d11B) & ~isnan(data.boron.d18O)),:);
data.input_parameters = jsondecode(fileread(data_directory+"/pH_Change/Input.json"));
data.alkalinity_constraints = jsondecode(fileread(data_directory+"/pH_Change/Alkalinity_Constraints.json"));
interpolation_ages = jsondecode(fileread(data_directory+"/Age/Interpolation_Age.json")).interpolation_ages;

%% Alkalinity
for alkalinity_shape_index = 1:size(data.alkalinity_constraints.shape,1)
    round_1.alkalinity.shape_samplers(alkalinity_shape_index) = Geochemistry_Helpers.Sampler(-1:0.0001:1,"Gaussian",data.alkalinity_constraints.shape(alkalinity_shape_index,2:3),'latin_hypercube_random').normalise();
    round_1.alkalinity.shape_samplers(alkalinity_shape_index).location = data.alkalinity_constraints.shape(alkalinity_shape_index,1);
end
round_1.alkalinity.scaling_sampler = Geochemistry_Helpers.Sampler(-12000:10:12000,"Flat",data.alkalinity_constraints.scale,"latin_hypercube_random").normalise();
round_1.alkalinity.scaling_sampler.getSamples(round_1.number_of_samples).shuffle();

round_1.alkalinity_gaussian_process = Geochemistry_Helpers.GaussianProcess("rbf",interpolation_ages);
round_1.alkalinity_gaussian_process.observations = round_1.alkalinity.shape_samplers;

round_1.alkalinity_gaussian_process.runKernel([0.5,0.05]);
round_1.alkalinity_gaussian_process.getSamples(round_1.number_of_samples);

%% Distributions and Samplers
% Generate samplers for each variable to represent uncertainties
% Initial state
round_1.co2.initial_sampler = Geochemistry_Helpers.Sampler(0:100e-6:10000e-6,"Flat",data.input_parameters.co2/1e6,"latin_hypercube_random").normalise();
round_1.saturation_state.initial_sampler = Geochemistry_Helpers.Sampler(0:0.1:12,"Flat",data.input_parameters.saturation_state,"latin_hypercube_random").normalise();

round_1.temperature.initial_samplers = [Geochemistry_Helpers.Sampler(-10:0.1:30,"Gaussian",data.input_parameters.initial_temperature(1,2:3),"latin_hypercube_random").normalise(data.input_parameters.initial_temperature(1,1)),Geochemistry_Helpers.Sampler(-10:0.1:30,"Gaussian",data.input_parameters.initial_temperature(2,2:3),"latin_hypercube_random").normalise(data.input_parameters.initial_temperature(2,1))];
round_1.temperature.initial_sampler = Geochemistry_Helpers.Sampler(-10:0.1:30,"Manual",round_1.temperature.initial_samplers(1).probabilities+round_1.temperature.initial_samplers(2).probabilities,"latin_hypercube_random").normalise();

% Constant across interval
round_1.calcium.sampler = Geochemistry_Helpers.Sampler(0:0.1:20,"Flat",data.input_parameters.calcium,"latin_hypercube_random").normalise();
round_1.magnesium.sampler = Geochemistry_Helpers.Sampler(20:0.1:61,"Flat",data.input_parameters.magnesium,"latin_hypercube_random").normalise();
round_1.species_calibration.gradient_sampler = Geochemistry_Helpers.Sampler(0:0.01:1.1,"Flat",data.input_parameters.species_calibration_gradient,"latin_hypercube_random").normalise();
round_1.species_calibration.intercept_sampler = Geochemistry_Helpers.Sampler(-20:0.1:20,"Flat",data.input_parameters.species_calibration_intercept,"latin_hypercube_random").normalise();

round_1.epsilon.samplers = [Geochemistry_Helpers.Sampler(23:0.01:30,"Gaussian",data.input_parameters.epsilon(1,2:3),'latin_hypercube_random').normalise(data.input_parameters.epsilon(1,1)),Geochemistry_Helpers.Sampler(23:0.01:30,"Gaussian",data.input_parameters.epsilon(2,2:3),'latin_hypercube_random').normalise(data.input_parameters.epsilon(2,1))];
round_1.epsilon.sampler = Geochemistry_Helpers.Sampler(23:0.01:30,"Manual",round_1.epsilon.samplers(1).probabilities+round_1.epsilon.samplers(2).probabilities,'latin_hypercube_random').normalise();

% Constant for all runs
round_1.oceanic_pressure = data.input_parameters.oceanic_pressure;
round_1.atmospheric_pressure = data.input_parameters.atmospheric_pressure;
round_1.salinity = data.input_parameters.salinity;

% Get the samples
round_1.co2.initial_sampler.getSamples(round_1.number_of_samples).shuffle();
round_1.saturation_state.initial_sampler.getSamples(round_1.number_of_samples).shuffle();
round_1.temperature.initial_sampler.getSamples(round_1.number_of_samples).shuffle();
round_1.calcium.sampler.getSamples(round_1.number_of_samples).shuffle();
round_1.magnesium.sampler.getSamples(round_1.number_of_samples).shuffle();

round_1.species_calibration.gradient_sampler.getSamples(round_1.number_of_samples).shuffle();
round_1.species_calibration.intercept_sampler.getSamples(round_1.number_of_samples).shuffle();
round_1.epsilon.sampler.getSamples(round_1.number_of_samples).shuffle();

% Time series
% d11B measured
round_1.d11B_measured.samplers = Geochemistry_Helpers.Sampler().create([height(data.boron),1]);
for d11B_sampler_index = 1:numel(round_1.d11B_measured.samplers)
    round_1.d11B_measured.samplers(d11B_sampler_index) = Geochemistry_Helpers.Sampler(-20:0.01:40,"Gaussian",[data.boron.d11B(d11B_sampler_index),data.boron.d11B_uncertainty(d11B_sampler_index)/2],'latin_hypercube_random').normalise();
    round_1.d11B_measured.samplers(d11B_sampler_index).location = data.boron.age(d11B_sampler_index);
end

round_1.d11B_measured.gaussian_process = Geochemistry_Helpers.GaussianProcess("rbf",interpolation_ages);
round_1.d11B_measured.gaussian_process.observations = round_1.d11B_measured.samplers';
round_1.d11B_measured.gaussian_process.runKernel([0.8,18000/1e6],linspace(1,5,numel(round_1.d11B_measured.samplers)));
round_1.d11B_measured.gaussian_process.getSamples(round_1.number_of_samples);

round_1.d11B_measured.samples = round_1.d11B_measured.gaussian_process.samples';

% d11B4
round_1.d11B_4.evolutions = BuCC.BoronSpeciesCalibration().create([numel(interpolation_ages),round_1.number_of_samples]);
for evolution_index = 1:size(round_1.d11B_4.evolutions,2)    
    % Species calibration
    round_1.d11B_4.evolutions(:,evolution_index).assignToAll("form","polynomial");
    round_1.d11B_4.evolutions(:,evolution_index).assignToAll("coefficients",[round_1.species_calibration.gradient_sampler.samples(evolution_index);round_1.species_calibration.intercept_sampler.samples(evolution_index)]');
    round_1.d11B_4.evolutions(:,evolution_index).d11B_measured.assignToEach("value",round_1.d11B_measured.samples(:,evolution_index));
    
    % Do species calculations
    round_1.d11B_4.evolutions(:,evolution_index).calculate();
end
round_1.d11B_4.samples = round_1.d11B_4.evolutions.d11B_4.value;

% Temperature
round_1.temperature.change_samplers = Geochemistry_Helpers.Sampler().create([height(data.boron.delta_temperature),1]);
for temperature_change_sampler_index = 1:numel(round_1.temperature.change_samplers)
    round_1.temperature.change_samplers(temperature_change_sampler_index) = Geochemistry_Helpers.Sampler(-20:0.1:60,"Gaussian",[data.boron.delta_temperature(temperature_change_sampler_index),data.boron.hansen_temperature_uncertainty(temperature_change_sampler_index)],'latin_hypercube_random').normalise();
    if temperature_change_sampler_index==1
        round_1.temperature.change_samplers(temperature_change_sampler_index) = Geochemistry_Helpers.Sampler(-20:0.1:20,"Gaussian",[data.boron.delta_temperature(temperature_change_sampler_index),0.01],'latin_hypercube_random').normalise();
    end
    round_1.temperature.change_samplers(temperature_change_sampler_index).location = data.boron.age(temperature_change_sampler_index);
end

round_1.temperature.change_gaussian_process = Geochemistry_Helpers.GaussianProcess("rbf",interpolation_ages);
round_1.temperature.change_gaussian_process.observations = round_1.temperature.change_samplers';
round_1.temperature.change_gaussian_process.runKernel([2,20000/1e6],1);
round_1.temperature.change_gaussian_process.getSamples(round_1.number_of_samples);
round_1.temperature.samples = round_1.temperature.initial_sampler.samples + round_1.temperature.change_gaussian_process.samples';

clearvars alkalinity_shape_index temperature_change_sampler_index evolution_index d11B_sampler_index

%% Initial pH
round_1.d11B_CO2.initial = BuCC.d11BCO2().create(round_1.number_of_samples);

% Species calibration
round_1.d11B_CO2.initial.species_calibration.assignToEach("coefficients",[round_1.species_calibration.gradient_sampler.samples;round_1.species_calibration.intercept_sampler.samples]');
round_1.d11B_CO2.initial.species_calibration.d11B_measured.assignToEach("value",round_1.d11B_measured.samples(1,:));

% Boron
round_1.d11B_CO2.initial.boron.pH.assignToAll("pValue",NaN);
round_1.d11B_CO2.initial.boron.d11B_sw.assignToAll("value",NaN);
round_1.d11B_CO2.initial.boron.assignToEach("epsilon",round_1.epsilon.sampler.samples);

% Carbonate chemistry
% Main parameters
round_1.d11B_CO2.initial.carbonate_chemistry.atmospheric_co2.assignToEach("partial_pressure",round_1.co2.initial_sampler.samples);
round_1.d11B_CO2.initial.carbonate_chemistry.assignToEach("saturation_state",round_1.saturation_state.initial_sampler.samples);
round_1.d11B_CO2.initial.carbonate_chemistry.assignToAll("units"," mol/kg");

% Ancillary
% Variable
round_1.d11B_CO2.initial.carbonate_chemistry.assignToEach("temperature",round_1.temperature.change_gaussian_process.samples(:,1));
round_1.d11B_CO2.initial.carbonate_chemistry.assignToEach("calcium",round_1.calcium.sampler.samples);
round_1.d11B_CO2.initial.carbonate_chemistry.assignToEach("magnesium",round_1.magnesium.sampler.samples);

% Constant
round_1.d11B_CO2.initial.carbonate_chemistry.assignToAll("salinity",round_1.salinity);
round_1.d11B_CO2.initial.carbonate_chemistry.assignToAll("oceanic_pressure",round_1.oceanic_pressure);
round_1.d11B_CO2.initial.carbonate_chemistry.assignToAll("atmospheric_pressure",round_1.atmospheric_pressure);

% Create a MyAMI object
myami = MyAMI.MyAMI("Precalculated",true);
round_1.d11B_CO2.initial.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

% Do carbonate chemistry calculations
round_1.d11B_CO2.initial.calculate();

% Collate the results of those calculations back to normal arrays
round_1.pH.initial_samples = round_1.d11B_CO2.initial.carbonate_chemistry.pH.pValue';
round_1.alkalinity.initial = round_1.d11B_CO2.initial.carbonate_chemistry.alkalinity;

round_1.alkalinity.samples = (((round_1.alkalinity_gaussian_process.samples.*repmat(round_1.alkalinity.scaling_sampler.samples',1,100))+round_1.alkalinity.initial*1e6)/1e6)';

%% We have starting conditions but can refine d11B_sw
% Collate the viable d11B_sw's and make a distribution
round_1.d11B_sw.samples = round_1.d11B_CO2.initial.boron.d11B_sw.value;
round_1.d11B_sw.sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.1:60,round_1.d11B_sw.samples,"latin_hypercube_random");

% We also have a d11B_sw constraint from the range of measured d11B
% So get the minimum and maximum
round_1.d11B_sw.minimum = max(round_1.d11B_4.samples);
round_1.d11B_sw.maximum_ratio = min(round_1.d11B_4.evolutions.d11B_4.ratio).*(1+round_1.epsilon.sampler.samples/1000);

round_1.d11B_sw.maximum_delta = Geochemistry_Helpers.delta().create(numel(round_1.d11B_sw.maximum_ratio));
round_1.d11B_sw.maximum_delta.assignToAll("standard","B");
round_1.d11B_sw.maximum_delta.assignToEach("ratio",round_1.d11B_sw.maximum_ratio);

round_1.d11B_sw.maximum = round_1.d11B_sw.maximum_delta.value';

%%
round_1.d11B_sw.valid = double(round_1.d11B_sw.samples>=round_1.d11B_sw.minimum' &  round_1.d11B_sw.samples<=round_1.d11B_sw.maximum');

round_2.number_of_samples = round_1.number_of_samples;
[round_2.d11B_sw.samples,round_2.resample_indices] = round_1.d11B_sw.sampler.resample(round_1.d11B_sw.valid,round_2.number_of_samples);

round_2.d11B_measured.samples = round_1.d11B_measured.samples(:,round_2.resample_indices);
round_2.d11B_4.samples = round_1.d11B_4.samples(:,round_2.resample_indices);
round_2.co2.initial_samples = round_1.co2.initial_sampler.samples(round_2.resample_indices);
round_2.saturation_state.initial_samples = round_1.saturation_state.initial_sampler.samples(round_2.resample_indices);
round_2.pH.initial_samples = round_1.pH.initial_samples(round_2.resample_indices);
round_2.calcium.samples = round_1.calcium.sampler.samples(round_2.resample_indices);
round_2.magnesium.samples = round_1.magnesium.sampler.samples(round_2.resample_indices);
round_2.epsilon.samples = round_1.epsilon.sampler.samples(round_2.resample_indices);
round_2.species_calibration.gradient_samples = round_1.species_calibration.gradient_sampler.samples(round_2.resample_indices);
round_2.species_calibration.intercept_samples = round_1.species_calibration.intercept_sampler.samples(round_2.resample_indices);
round_2.temperature.samples = round_1.temperature.samples(:,round_2.resample_indices);
round_2.alkalinity_samples = round_1.alkalinity.initial(round_2.resample_indices);

round_2.d11B_sw.sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.1:60,round_2.d11B_sw.samples,'latin_hypercube');
round_2.co2.initial_sampler = Geochemistry_Helpers.Sampler.fromSamples(0:100e-6:10000e-6,round_2.co2.initial_samples,'latin_hypercube');
round_2.saturation_state.initial_sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.1:12,round_2.saturation_state.initial_samples,'latin_hypercube');
round_2.calcium.sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.1:20,round_2.calcium.samples,'latin_hypercube');
round_2.magnesium.sampler = Geochemistry_Helpers.Sampler.fromSamples(20:0.1:61,round_2.magnesium.samples,'latin_hypercube');
round_2.epsilon.sampler = Geochemistry_Helpers.Sampler.fromSamples(23:0.01:30,round_2.epsilon.samples,'latin_hypercube');
round_2.species_calibration.gradient_sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.01:1.1,round_2.species_calibration.gradient_samples,'latin_hypercube');
round_2.species_calibration.intercept_sampler = Geochemistry_Helpers.Sampler.fromSamples(-10:0.1:10,round_2.species_calibration.intercept_samples,'latin_hypercube');

round_2.alkalinity.samples = round_1.alkalinity.samples;

round_2.salinity = round_1.salinity;
round_2.oceanic_pressure = round_1.oceanic_pressure;
round_2.atmospheric_pressure = round_1.atmospheric_pressure;

%%
round_2.d11B_CO2.evolutions = BuCC.d11BCO2().create([numel(interpolation_ages),round_2.number_of_samples,2]);
for evolution_index = 1:size(round_2.d11B_CO2.evolutions,2)
    % Low omega
    % Species calibration
    round_2.d11B_CO2.evolutions(:,evolution_index,1).species_calibration.assignToAll("coefficients",[round_2.species_calibration.gradient_samples(evolution_index);round_2.species_calibration.intercept_samples(evolution_index)]');
    round_2.d11B_CO2.evolutions(:,evolution_index,1).species_calibration.d11B_measured.assignToEach("value",round_2.d11B_measured.samples(:,evolution_index));
    
    % Boron
    round_2.d11B_CO2.evolutions(:,evolution_index,1).boron.pH.assignToAll("pValue",NaN);
    round_2.d11B_CO2.evolutions(:,evolution_index,1).boron.d11B_sw.assignToAll("value",round_2.d11B_sw.samples(evolution_index));
    round_2.d11B_CO2.evolutions(:,evolution_index,1).boron.assignToAll("epsilon",round_2.epsilon.samples(evolution_index));
    
    % Carbonate chemistry
    % Main parameters
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",NaN);
    
    round_2.saturation_state.lower_bound = 1*ones(numel(interpolation_ages),1);
    round_2.saturation_state.lower_bound(interpolation_ages>data.boron.age(9)) = 5;
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToEach("saturation_state",round_2.saturation_state.lower_bound);
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("units"," mol/kg");
    
    % Ancillary
    % Variable
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToEach("temperature",round_2.temperature.samples(:,evolution_index));
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("calcium",round_2.calcium.samples(evolution_index));
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("magnesium",round_2.magnesium.samples(evolution_index));
    
    % Constant
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("salinity",round_2.salinity);
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("oceanic_pressure",round_2.oceanic_pressure);
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("atmospheric_pressure",round_2.atmospheric_pressure);
    
    % Create a MyAMI object
%     myami = MyAMI.MyAMI("Precalculated",true);
    round_2.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);
    
    % Do carbonate chemistry calculations
    round_2.d11B_CO2.evolutions(:,evolution_index,1).calculate();
    
    % High omega
    % Species calibration
    round_2.d11B_CO2.evolutions(:,evolution_index,2).species_calibration.assignToAll("coefficients",[round_2.species_calibration.gradient_samples(evolution_index);round_2.species_calibration.intercept_samples(evolution_index)]');
    round_2.d11B_CO2.evolutions(:,evolution_index,2).species_calibration.d11B_measured.assignToEach("value",round_2.d11B_measured.samples(:,evolution_index));
    
    % Boron
    round_2.d11B_CO2.evolutions(:,evolution_index,2).boron.pH.assignToAll("pValue",NaN);
    round_2.d11B_CO2.evolutions(:,evolution_index,2).boron.d11B_sw.assignToAll("value",round_2.d11B_sw.samples(evolution_index));
    round_2.d11B_CO2.evolutions(:,evolution_index,2).boron.assignToAll("epsilon",round_2.epsilon.samples(evolution_index));
    
    % Carbonate chemistry
    % Main parameters
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",NaN);
    round_2.saturation_state.upper_bound = 12;
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("saturation_state",round_2.saturation_state.upper_bound);
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("units"," mol/kg");
    
    % Ancillary
    % Variable
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.assignToEach("temperature",round_2.temperature.samples(:,evolution_index));
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("calcium",round_2.calcium.samples(evolution_index));
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("magnesium",round_2.magnesium.samples(evolution_index));
    
    % Constant
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("salinity",round_2.salinity);
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("oceanic_pressure",round_2.oceanic_pressure);
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.assignToAll("atmospheric_pressure",round_2.atmospheric_pressure);
    
    % Create a MyAMI object
%     myami = MyAMI.MyAMI("Precalculated",true);
    round_2.d11B_CO2.evolutions(:,evolution_index,2).carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);
    
    % Do carbonate chemistry calculations
    round_2.d11B_CO2.evolutions(:,evolution_index,2).calculate();  
end

% Collate the results of those calculations back to normal arrays
round_2.pH.value = round_2.d11B_CO2.evolutions.carbonate_chemistry.pH.pValue;
round_2.alkalinity.boundaries = round_2.d11B_CO2.evolutions.carbonate_chemistry.alkalinity;

%% To calculate evolutions we need to know how alkalinity may have changed
% Do this using Gaussian process

round_2.alkalinity.viable = all(round_2.alkalinity.samples>=round_2.alkalinity.boundaries(:,:,1)) & all(round_2.alkalinity.samples<=round_2.alkalinity.boundaries(:,:,2));

round_3.number_of_samples = sum(round_2.alkalinity.viable);

round_3.co2.initial_samples = round_2.co2.initial_samples(:,round_2.alkalinity.viable);
round_3.saturation_state.initial_samples = round_2.saturation_state.initial_samples(:,round_2.alkalinity.viable);
round_3.alkalinity.samples = round_2.alkalinity.samples(:,round_2.alkalinity.viable);
round_3.pH.initial_samples = round_2.pH.initial_samples(:,round_2.alkalinity.viable);
round_3.pH.samples = round_2.d11B_CO2.evolutions(:,:,1).carbonate_chemistry.pH.pValue(:,round_2.alkalinity.viable);
round_3.d11B_measured.samples = round_2.d11B_measured.samples(:,round_2.alkalinity.viable);
round_3.d11B_4.samples = round_2.d11B_4.samples(:,round_2.resample_indices);
round_3.species_calibration.gradient_samples = round_2.species_calibration.gradient_samples(round_2.alkalinity.viable);
round_3.species_calibration.intercept_samples = round_2.species_calibration.intercept_samples(round_2.alkalinity.viable);

round_3.d11B_sw.samples = round_2.d11B_sw.samples(round_2.alkalinity.viable);
round_3.epsilon.samples = round_2.epsilon.samples(round_2.alkalinity.viable);
round_3.calcium.samples = round_2.calcium.samples(round_2.alkalinity.viable);
round_3.magnesium.samples = round_2.magnesium.samples(round_2.alkalinity.viable);

round_3.temperature.samples = round_2.temperature.samples(:,round_2.alkalinity.viable);

round_3.d11B_sw.sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.1:60,round_3.d11B_sw.samples,'latin_hypercube');
round_3.co2.initial_sampler = Geochemistry_Helpers.Sampler.fromSamples(0:100e-6:10000e-6,round_3.co2.initial_samples,'latin_hypercube');
round_3.saturation_state.initial_sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.1:12,round_3.saturation_state.initial_samples,'latin_hypercube');
round_3.calcium.sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.1:20,round_3.calcium.samples,'latin_hypercube');
round_3.magnesium.sampler = Geochemistry_Helpers.Sampler.fromSamples(20:0.1:61,round_3.magnesium.samples,'latin_hypercube');
round_3.epsilon.sampler = Geochemistry_Helpers.Sampler.fromSamples(23:0.01:30,round_3.epsilon.samples,'latin_hypercube');
round_3.species_calibration.gradient_sampler = Geochemistry_Helpers.Sampler.fromSamples(0:0.01:1.1,round_3.species_calibration.gradient_samples,'latin_hypercube');
round_3.species_calibration.intercept_sampler = Geochemistry_Helpers.Sampler.fromSamples(-10:0.1:10,round_3.species_calibration.intercept_samples,'latin_hypercube');

round_3.salinity = round_2.salinity;
round_3.oceanic_pressure = round_2.oceanic_pressure;
round_3.atmospheric_pressure = round_2.atmospheric_pressure;

%%
round_3.d11B_CO2.evolutions = BuCC.d11BCO2().create([numel(interpolation_ages),size(round_3.alkalinity.samples,2)]);
for evolution_index = 1:size(round_3.alkalinity.samples,2)
    % Species calibration
    round_3.d11B_CO2.evolutions(:,evolution_index,1).species_calibration.assignToAll("coefficients",[round_3.species_calibration.gradient_samples(evolution_index),round_3.species_calibration.intercept_samples(evolution_index)]');
    round_3.d11B_CO2.evolutions(:,evolution_index,1).species_calibration.d11B_measured.assignToEach("value",round_3.d11B_measured.samples(:,evolution_index));
    
    % Boron
    round_3.d11B_CO2.evolutions(:,evolution_index,1).boron.pH.assignToAll("pValue",NaN);
    round_3.d11B_CO2.evolutions(:,evolution_index,1).boron.d11B_sw.assignToAll("value",round_3.d11B_sw.samples(evolution_index));
    round_3.d11B_CO2.evolutions(:,evolution_index,1).boron.assignToAll("epsilon",round_3.epsilon.samples(evolution_index));
    
    % Carbonate chemistry
    % Main parameters
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",NaN);
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToEach("alkalinity",round_3.alkalinity.samples(:,evolution_index));
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("units"," mol/kg");
    
    % Ancillary
    % Variable
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToEach("temperature",round_3.temperature.samples(:,evolution_index));
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("calcium",round_3.calcium.samples(evolution_index));
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("magnesium",round_3.magnesium.samples(evolution_index));
    
    % Constant
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("salinity",round_3.salinity);
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("oceanic_pressure",round_3.oceanic_pressure);
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.assignToAll("atmospheric_pressure",round_3.atmospheric_pressure);
    
    % Create a MyAMI object
%     myami = MyAMI.MyAMI("Precalculated",true);
    round_3.d11B_CO2.evolutions(:,evolution_index,1).carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);
    
    % Do carbonate chemistry calculations
    round_3.d11B_CO2.evolutions(:,evolution_index).calculate();
end

round_3.co2.samples = round_3.d11B_CO2.evolutions.carbonate_chemistry.atmospheric_co2.x;
round_3.saturation_state.samples = round_3.d11B_CO2.evolutions.carbonate_chemistry.saturation_state;
round_3.alkalinity.samples = round_3.d11B_CO2.evolutions.carbonate_chemistry.alkalinity;
round_3.pH.samples = round_3.d11B_CO2.evolutions.carbonate_chemistry.pH.pValue;
round_3.d11B_sw.samples = round_3.d11B_CO2.evolutions.boron.d11B_sw.value;
