% Minimum_pH_Change_Variation
% This script quantifies the change in pH with a variation in each of the
% main parameters needed: calcium, magnesium, temperature change, initial
% temperature, saturation state and CO2

% Each block below follows the same pattern:
% - Set up the variables, with one as an array of options
% - Use the preperturbation value for the initial pH with the various
% values for the parameter in the array
% - Use the postperturbation value for the minimum pH with the various
% values for the parameter in the array
% - Do the difference

% The results are then saved to a JSON file

clear

%% Load Data
data_directory = "./../../../Data/";
temperature_data = readtable(data_directory+"/Temperature/TJ_d18O_d13C.xlsx","Sheet","Averaged");
boron_data = readtable(data_directory+"/Boron/TJ_d11B.xlsx");

% Create a map to hold the results
pH_difference = containers.Map();

%% Default parameters
default.saturation_state = 10.7; % From Ridgwell modelling
default.co2 = 500e-6; % From Witkowski
default.initial_temperature = 25.9; %
default.temperature_change = 3.4;
default.salinity = 35; % Assumed
default.pressure = 0; % Assumed?
default.calcium = 10; % Horita
default.magnesium = 52; % Horita
default.epsilon = 27.8;

range.saturation_state = 10.4:0.1:11; % From Ridgwell modelling
range.co2 = 400e-6:50e-6:600e-6; % From Witkowski
range.initial_temperature = 20:0.1:31.8; % From bulk temperatures
range.temperature_change = 2.4:0.2:4.4;
range.calcium = 9:0.5:11; % Horita
range.magnesium = 51:0.5:53; % Horita
range.epsilon = 27.5:0.1:28.1;

%%
background.data = boron_data(1:9,:);
background.d11B4 = mean(background.data.d11B);

perturbation.data = boron_data(10:end,:);
perturbation.data_sorted = sortrows(perturbation.data,"d11B");
perturbation.d11B4 = mean(perturbation.data_sorted.d11B(1:2));

d11B4_change = background.d11B4-perturbation.d11B4;

%% Ca
% Process background
background.samples = BuCC.d11BCO2().create(numel(range.calcium));
background.samples.collate("species_calibration").collate("d11B_measured").assignToAll("value",background.d11B4);

background.samples.boron.pH.assignToAll("pValue",NaN);
background.samples.boron.assignToAll("epsilon",default.epsilon);
background.samples.boron.d11B_sw.assignToAll("value",NaN);

background.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature);
background.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
background.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
background.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
background.samples.carbonate_chemistry.assignToEach("calcium",range.calcium);
background.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

background.samples.carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",default.co2);
background.samples.carbonate_chemistry.assignToAll("saturation_state",default.saturation_state);

background.samples.calculate();

% Process perturbation
perturbation.samples = BuCC.d11BCO2().create(numel(range.calcium));
perturbation.samples.species_calibration.d11B_measured.assignToAll("value",perturbation.d11B4);

perturbation.samples.boron.pH.assignToAll("pValue",NaN);
perturbation.samples.boron.assignToAll("epsilon",default.epsilon);
perturbation.samples.boron.d11B_sw.assignToEach("value",background.samples.collate("boron").collate("d11B_sw").collate("value"));

perturbation.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature+default.temperature_change);
perturbation.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
perturbation.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
perturbation.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
perturbation.samples.carbonate_chemistry.assignToEach("calcium",range.calcium);
perturbation.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

perturbation.samples.carbonate_chemistry.assignToEach("alkalinity",background.samples.collate("carbonate_chemistry").collate("alkalinity"));

perturbation.samples.calculate();

% Calculate change
background.pH = background.samples.carbonate_chemistry.pH.pValue;
perturbation.pH = perturbation.samples.carbonate_chemistry.pH.pValue;

pH_change = perturbation.pH-background.pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("calcium") = pH_change;

%% Magnesium
% Process background
background.samples = BuCC.d11BCO2().create(numel(range.magnesium));
background.samples.species_calibration.d11B_measured.assignToAll("value",background.d11B4);

background.samples.boron.pH.assignToAll("pValue",NaN);
background.samples.boron.assignToAll("epsilon",default.epsilon);
background.samples.boron.d11B_sw.assignToAll("value",NaN);

background.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature);
background.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
background.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
background.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
background.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
background.samples.carbonate_chemistry.assignToEach("magnesium",range.magnesium);

background.samples.carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",default.co2);
background.samples.carbonate_chemistry.assignToAll("saturation_state",default.saturation_state);

background.samples.calculate();

% Process perturbation
perturbation.samples = BuCC.d11BCO2().create(numel(range.magnesium));
perturbation.samples.species_calibration.d11B_measured.assignToAll("value",perturbation.d11B4);

perturbation.samples.boron.pH.assignToAll("pValue",NaN);
perturbation.samples.boron.assignToAll("epsilon",default.epsilon);
perturbation.samples.boron.d11B_sw.assignToEach("value",background.samples.collate("boron").collate("d11B_sw").collate("value"));

perturbation.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature+default.temperature_change);
perturbation.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
perturbation.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
perturbation.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
perturbation.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
perturbation.samples.carbonate_chemistry.assignToEach("magnesium",range.magnesium);

perturbation.samples.carbonate_chemistry.assignToEach("alkalinity",background.samples.collate("carbonate_chemistry").collate("alkalinity"));

perturbation.samples.calculate();

% Calculate change
background.pH = background.samples.carbonate_chemistry.pH.pValue;
perturbation.pH = perturbation.samples.carbonate_chemistry.pH.pValue;

pH_change = perturbation.pH-background.pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("magnesium") = pH_change;

%% Temperature change
% Process background
background.samples = BuCC.d11BCO2().create(numel(range.temperature_change));
background.samples.species_calibration.d11B_measured.assignToAll("value",background.d11B4);

background.samples.boron.pH.assignToAll("pValue",NaN);
background.samples.boron.assignToAll("epsilon",default.epsilon);
background.samples.boron.d11B_sw.assignToAll("value",NaN);

background.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature);
background.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
background.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
background.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
background.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
background.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

background.samples.carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",default.co2);
background.samples.carbonate_chemistry.assignToAll("saturation_state",default.saturation_state);

background.samples.calculate();

% Process perturbation
perturbation.samples = BuCC.d11BCO2().create(numel(range.temperature_change));
perturbation.samples.species_calibration.d11B_measured.assignToAll("value",perturbation.d11B4);

perturbation.samples.boron.pH.assignToAll("pValue",NaN);
perturbation.samples.boron.assignToAll("epsilon",default.epsilon);
perturbation.samples.boron.d11B_sw.assignToEach("value",background.samples.boron.d11B_sw.value);

perturbation.samples.carbonate_chemistry.assignToEach("temperature",default.initial_temperature+range.temperature_change);
perturbation.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
perturbation.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
perturbation.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
perturbation.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
perturbation.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

perturbation.samples.carbonate_chemistry.assignToEach("alkalinity",background.samples.carbonate_chemistry.alkalinity);

perturbation.samples.calculate();

% Calculate change
background.pH = background.samples.carbonate_chemistry.pH.pValue;
perturbation.pH = perturbation.samples.carbonate_chemistry.pH.pValue;

pH_change = perturbation.pH-background.pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("temperature_change") = pH_change;

%% Temperature start
% Process background
background.samples = BuCC.d11BCO2().create(numel(range.initial_temperature));
background.samples.species_calibration.d11B_measured.assignToAll("value",background.d11B4);

background.samples.boron.pH.assignToAll("pValue",NaN);
background.samples.boron.assignToAll("epsilon",default.epsilon);
background.samples.boron.d11B_sw.assignToAll("value",NaN);

background.samples.carbonate_chemistry.assignToEach("temperature",range.initial_temperature);
background.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
background.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
background.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
background.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
background.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

background.samples.carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",default.co2);
background.samples.carbonate_chemistry.assignToAll("saturation_state",default.saturation_state);

background.samples.calculate();

%
perturbation.samples = BuCC.d11BCO2().create(numel(range.initial_temperature));
perturbation.samples.species_calibration.d11B_measured.assignToAll("value",perturbation.d11B4);

perturbation.samples.boron.pH.assignToAll("pValue",NaN);
perturbation.samples.boron.assignToAll("epsilon",default.epsilon);
perturbation.samples.boron.d11B_sw.assignToEach("value",background.samples.boron.d11B_sw.value);

perturbation.samples.carbonate_chemistry.assignToEach("temperature",range.initial_temperature+default.temperature_change);
perturbation.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
perturbation.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
perturbation.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
perturbation.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
perturbation.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

perturbation.samples.carbonate_chemistry.assignToEach("alkalinity",background.samples.carbonate_chemistry.alkalinity);

perturbation.samples.calculate();

% Calculate change
background.pH = background.samples.carbonate_chemistry.pH.pValue;
perturbation.pH = perturbation.samples.carbonate_chemistry.pH.pValue;

pH_change = perturbation.pH-background.pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("initial_temperature") = pH_change;

%% Saturation State
% Process background
background.samples = BuCC.d11BCO2().create(numel(range.saturation_state));
background.samples.species_calibration.d11B_measured.assignToAll("value",background.d11B4);

background.samples.boron.pH.assignToAll("pValue",NaN);
background.samples.boron.assignToAll("epsilon",default.epsilon);
background.samples.boron.d11B_sw.assignToAll("value",NaN);

background.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature);
background.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
background.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
background.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
background.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
background.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

background.samples.carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",default.co2);
background.samples.carbonate_chemistry.assignToEach("saturation_state",range.saturation_state);

background.samples.calculate();

%
perturbation.samples = BuCC.d11BCO2().create(numel(range.saturation_state));
perturbation.samples.species_calibration.d11B_measured.assignToAll("value",perturbation.d11B4);

perturbation.samples.boron.pH.assignToAll("pValue",NaN);
perturbation.samples.boron.assignToAll("epsilon",default.epsilon);
perturbation.samples.boron.d11B_sw.assignToEach("value",background.samples.boron.d11B_sw.value);

perturbation.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature+default.temperature_change);
perturbation.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
perturbation.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
perturbation.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
perturbation.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
perturbation.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

perturbation.samples.carbonate_chemistry.assignToEach("alkalinity",background.samples.carbonate_chemistry.alkalinity);

perturbation.samples.calculate();

% Calculate change
background.pH = background.samples.carbonate_chemistry.pH.pValue;
perturbation.pH = perturbation.samples.carbonate_chemistry.pH.pValue;

pH_change = perturbation.pH-background.pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("saturation_state") = pH_change;

%% CO2
% Process background
background.samples = BuCC.d11BCO2().create(numel(range.co2));
background.samples.species_calibration.d11B_measured.assignToAll("value",background.d11B4);

background.samples.boron.pH.assignToAll("pValue",NaN);
background.samples.boron.assignToAll("epsilon",default.epsilon);
background.samples.boron.d11B_sw.assignToAll("value",NaN);

background.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature);
background.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
background.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
background.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
background.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
background.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

background.samples.carbonate_chemistry.atmospheric_co2.assignToEach("partial_pressure",range.co2);
background.samples.carbonate_chemistry.assignToAll("saturation_state",default.saturation_state);

background.samples.calculate();

% Process perturbation
perturbation.samples = BuCC.d11BCO2().create(numel(range.co2));
perturbation.samples.species_calibration.d11B_measured.assignToAll("value",perturbation.d11B4);

perturbation.samples.boron.pH.assignToAll("pValue",NaN);
perturbation.samples.boron.assignToAll("epsilon",default.epsilon);
perturbation.samples.boron.d11B_sw.assignToEach("value",background.samples.boron.d11B_sw.value);

perturbation.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature+default.temperature_change);
perturbation.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
perturbation.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
perturbation.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
perturbation.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
perturbation.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

perturbation.samples.carbonate_chemistry.assignToEach("alkalinity",background.samples.carbonate_chemistry.alkalinity);

perturbation.samples.calculate();

% Calculate change
background.pH = background.samples.carbonate_chemistry.pH.pValue;
perturbation.pH = perturbation.samples.carbonate_chemistry.pH.pValue;

pH_change = perturbation.pH-background.pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("co2") = pH_change;

%% Epsilon
% Process background
background.samples = BuCC.d11BCO2().create(numel(range.epsilon));
background.samples.species_calibration.d11B_measured.assignToAll("value",background.d11B4);

background.samples.boron.pH.assignToAll("pValue",NaN);
background.samples.boron.assignToEach("epsilon",range.epsilon);
background.samples.boron.d11B_sw.assignToAll("value",NaN);

background.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature);
background.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
background.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
background.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
background.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
background.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

background.samples.carbonate_chemistry.atmospheric_co2.assignToAll("partial_pressure",default.co2);
background.samples.carbonate_chemistry.assignToAll("saturation_state",default.saturation_state);

background.samples.calculate();

% Process perturbation
perturbation.samples = BuCC.d11BCO2().create(numel(range.epsilon));
perturbation.samples.species_calibration.d11B_measured.assignToAll("value",perturbation.d11B4);

perturbation.samples.boron.pH.assignToAll("pValue",NaN);
perturbation.samples.boron.assignToEach("epsilon",range.epsilon);
perturbation.samples.boron.d11B_sw.assignToEach("value",background.samples.boron.d11B_sw.value);

perturbation.samples.carbonate_chemistry.assignToAll("temperature",default.initial_temperature+default.temperature_change);
perturbation.samples.carbonate_chemistry.assignToAll("salinity",default.salinity);
perturbation.samples.carbonate_chemistry.assignToAll("oceanic_pressure",default.pressure);
perturbation.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
perturbation.samples.carbonate_chemistry.assignToAll("calcium",default.calcium);
perturbation.samples.carbonate_chemistry.assignToAll("magnesium",default.magnesium);

perturbation.samples.carbonate_chemistry.assignToEach("alkalinity",background.samples.carbonate_chemistry.alkalinity);

perturbation.samples.calculate();

% Calculate change
background.pH = background.samples.carbonate_chemistry.pH.pValue;
perturbation.pH = perturbation.samples.carbonate_chemistry.pH.pValue;

pH_change = perturbation.pH-background.pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("epsilon") = pH_change;

%% Saving
output_file = fopen(data_directory+"/Minimum_pH_Change/Individual_Sensitivity.json","w");

to_encode = [];
to_encode.calcium = range.calcium;
to_encode.pH_difference = pH_difference("calcium");
json = jsonencode(to_encode);
fprintf(output_file,"["+json+",");

to_encode = [];
to_encode.magnesium = range.magnesium;
to_encode.pH_difference = pH_difference("magnesium");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.temperature_change = range.temperature_change;
to_encode.pH_difference = pH_difference("temperature_change");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.initial_temperature = range.initial_temperature;
to_encode.pH_difference = pH_difference("initial_temperature");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.saturation_state = range.saturation_state;
to_encode.pH_difference = pH_difference("saturation_state");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.co2 = range.co2;
to_encode.pH_difference = pH_difference("co2");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.epsilon = range.epsilon;
to_encode.pH_difference = pH_difference("epsilon");
json = jsonencode(to_encode);
fprintf(output_file,json+"]");

fclose(output_file);

