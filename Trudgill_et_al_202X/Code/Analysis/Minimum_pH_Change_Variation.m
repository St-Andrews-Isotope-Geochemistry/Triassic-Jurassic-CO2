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

%% Load Data
temperature_data = readtable("./../../Data/TJ_Temperature.xlsx");
boron_data = readtable("./../../Data/TJ_d11B_pH.xlsx");

% Create a map to hold the results
pH_difference = containers.Map();

%% Ca
saturation_state_maximum = 10.7; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 21.4; % From bulk temperatures
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = 17; % Horita
mg = 28; % Horita

cas = [10:1:24];

%
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4_minimum = nanmin(perturbation_data.d11B);

d11B4_change = background_d11B4-perturbation_d11B4_minimum;

% Process with original bulk temperature for comparison
initial = BuCC.d11BCO2().create(numel(cas));
initial.collate("species_calibration").collate("d11B_measured").assignToAll("value",background_d11B4);

initial.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial.collate("boron").assignToAll("epsilon",27.2);
initial.collate("boron").collate("d11B_sw").assignToAll("value",NaN);

cc = initial.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToEach("calcium",cas);
cc.assignToAll("magnesium",mg);

initial.collate("carbonate_chemistry").collate("atmospheric_co2").assignToAll("partial_pressure",co2_minimum);
initial.collate("carbonate_chemistry").assignToAll("saturation_state",saturation_state_maximum);

initial.calculate();

%
after = BuCC.d11BCO2().create(numel(cas));
after.collate("species_calibration").collate("d11B_measured").assignToAll("value",perturbation_d11B4_minimum);

after.collate("boron").collate("pH").assignToAll("pValue",NaN);
after.collate("boron").assignToAll("epsilon",27.2);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToEach("calcium",cas);
cc.assignToAll("magnesium",mg);

after.collate("carbonate_chemistry").assignToEach("alkalinity",initial.collate("carbonate_chemistry").collate("alkalinity"));

after.calculate();

%
before_pH = initial.collate("carbonate_chemistry").collate("pH").collate("pValue");
after_pH = after.collate("carbonate_chemistry").collate("pH").collate("pValue");

pH_difference("calcium") = after_pH-before_pH;

%% Magnesium
saturation_state_maximum = 10.7; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 21.4; % From bulk temperatures
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = 17; % Horita
mg = 28; % Horita

mgs = [20:1:36];

%
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4_minimum = nanmin(perturbation_data.d11B);

d11B4_change = background_d11B4-perturbation_d11B4_minimum;

% Process with original bulk temperature for comparison
initial = BuCC.d11BCO2().create(numel(mgs));
initial.collate("species_calibration").collate("d11B_measured").assignToAll("value",background_d11B4);

initial.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial.collate("boron").assignToAll("epsilon",27.2);
initial.collate("boron").collate("d11B_sw").assignToAll("value",NaN);

cc = initial.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToEach("magnesium",mgs);

initial.collate("carbonate_chemistry").collate("atmospheric_co2").assignToAll("partial_pressure",co2_minimum);
initial.collate("carbonate_chemistry").assignToAll("saturation_state",saturation_state_maximum);

initial.calculate();

%
after = BuCC.d11BCO2().create(numel(mgs));
after.collate("species_calibration").collate("d11B_measured").assignToAll("value",perturbation_d11B4_minimum);

after.collate("boron").collate("pH").assignToAll("pValue",NaN);
after.collate("boron").assignToAll("epsilon",27.2);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToEach("magnesium",mgs);

after.collate("carbonate_chemistry").assignToEach("alkalinity",initial.collate("carbonate_chemistry").collate("alkalinity"));

after.calculate();

%
before_pH = initial.collate("carbonate_chemistry").collate("pH").collate("pValue");
after_pH = after.collate("carbonate_chemistry").collate("pH").collate("pValue");

pH_difference("magnesium") = after_pH-before_pH;

%% Temperature change
saturation_state_maximum = 10.7; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 21.4; % From bulk temperatures
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = 17; % Horita
mg = 28; % Horita

initial_temperature = 20;
temperature_changes = [-20:1:20];

%
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4_minimum = nanmin(perturbation_data.d11B);

d11B4_change = background_d11B4-perturbation_d11B4_minimum;

% Process with original bulk temperature for comparison
initial = BuCC.d11BCO2().create(numel(temperature_changes));
initial.collate("species_calibration").collate("d11B_measured").assignToAll("value",background_d11B4);

initial.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial.collate("boron").assignToAll("epsilon",27.2);
initial.collate("boron").collate("d11B_sw").assignToAll("value",NaN);

cc = initial.collate("carbonate_chemistry");
cc.assignToAll("temperature",initial_temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

initial.collate("carbonate_chemistry").collate("atmospheric_co2").assignToAll("partial_pressure",co2_minimum);
initial.collate("carbonate_chemistry").assignToAll("saturation_state",saturation_state_maximum);

initial.calculate();

%
after = BuCC.d11BCO2().create(numel(temperature_changes));
after.collate("species_calibration").collate("d11B_measured").assignToAll("value",perturbation_d11B4_minimum);

after.collate("boron").collate("pH").assignToAll("pValue",NaN);
after.collate("boron").assignToAll("epsilon",27.2);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToEach("temperature",initial_temperature+temperature_changes);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

after.collate("carbonate_chemistry").assignToEach("alkalinity",initial.collate("carbonate_chemistry").collate("alkalinity"));

after.calculate();

%
before_pH = initial.collate("carbonate_chemistry").collate("pH").collate("pValue");
after_pH = after.collate("carbonate_chemistry").collate("pH").collate("pValue");

pH_difference("temperature_change") = after_pH-before_pH;

%% Temperature start
saturation_state_maximum = 10.7; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 21.4; % From bulk temperatures
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = 17; % Horita
mg = 28; % Horita

initial_temperatures = [1.4:1:41.4];

%
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4_minimum = nanmin(perturbation_data.d11B);

d11B4_change = background_d11B4-perturbation_d11B4_minimum;

% Process with original bulk temperature for comparison
initial = BuCC.d11BCO2().create(numel(initial_temperatures));
initial.collate("species_calibration").collate("d11B_measured").assignToAll("value",background_d11B4);

initial.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial.collate("boron").assignToAll("epsilon",27.2);
initial.collate("boron").collate("d11B_sw").assignToAll("value",NaN);

cc = initial.collate("carbonate_chemistry");
cc.assignToEach("temperature",initial_temperatures);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

initial.collate("carbonate_chemistry").collate("atmospheric_co2").assignToAll("partial_pressure",co2_minimum);
initial.collate("carbonate_chemistry").assignToAll("saturation_state",saturation_state_maximum);

initial.calculate();

%
after = BuCC.d11BCO2().create(numel(initial_temperatures));
after.collate("species_calibration").collate("d11B_measured").assignToAll("value",perturbation_d11B4_minimum);

after.collate("boron").collate("pH").assignToAll("pValue",NaN);
after.collate("boron").assignToAll("epsilon",27.2);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToEach("temperature",initial_temperatures);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

after.collate("carbonate_chemistry").assignToEach("alkalinity",initial.collate("carbonate_chemistry").collate("alkalinity"));

after.calculate();

%
before_pH = initial.collate("carbonate_chemistry").collate("pH").collate("pValue");
after_pH = after.collate("carbonate_chemistry").collate("pH").collate("pValue");

pH_difference("initial_temperature") = after_pH-before_pH;

%% Saturation State
saturation_state_maximum = 10.7; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 21.4; % From bulk temperatures
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = 17; % Horita
mg = 28; % Horita

saturation_states = [8.7:0.1:12.7];

%
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4_minimum = nanmin(perturbation_data.d11B);

d11B4_change = background_d11B4-perturbation_d11B4_minimum;

% Process with original bulk temperature for comparison
initial = BuCC.d11BCO2().create(numel(saturation_states));
initial.collate("species_calibration").collate("d11B_measured").assignToAll("value",background_d11B4);

initial.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial.collate("boron").assignToAll("epsilon",27.2);
initial.collate("boron").collate("d11B_sw").assignToAll("value",NaN);

cc = initial.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

initial.collate("carbonate_chemistry").collate("atmospheric_co2").assignToAll("partial_pressure",co2_minimum);
initial.collate("carbonate_chemistry").assignToEach("saturation_state",saturation_states);

initial.calculate();

%
after = BuCC.d11BCO2().create(numel(saturation_states));
after.collate("species_calibration").collate("d11B_measured").assignToAll("value",perturbation_d11B4_minimum);

after.collate("boron").collate("pH").assignToAll("pValue",NaN);
after.collate("boron").assignToAll("epsilon",27.2);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

after.collate("carbonate_chemistry").assignToEach("alkalinity",initial.collate("carbonate_chemistry").collate("alkalinity"));

after.calculate();

%
before_pH = initial.collate("carbonate_chemistry").collate("pH").collate("pValue");
after_pH = after.collate("carbonate_chemistry").collate("pH").collate("pValue");

pH_difference("saturation_state") = after_pH-before_pH;


%% CO2
saturation_state_maximum = 10.7; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 21.4; % From bulk temperatures
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = 17; % Horita
mg = 28; % Horita

co2s = [300e-6:10e-6:700e-6];

%
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4_minimum = nanmin(perturbation_data.d11B);

d11B4_change = background_d11B4-perturbation_d11B4_minimum;

% Process with original bulk temperature for comparison
initial = BuCC.d11BCO2().create(numel(co2s));
initial.collate("species_calibration").collate("d11B_measured").assignToAll("value",background_d11B4);

initial.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial.collate("boron").assignToAll("epsilon",27.2);
initial.collate("boron").collate("d11B_sw").assignToAll("value",NaN);

cc = initial.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

initial.collate("carbonate_chemistry").collate("atmospheric_co2").assignToEach("partial_pressure",co2s);
initial.collate("carbonate_chemistry").assignToAll("saturation_state",saturation_state_maximum);

initial.calculate();

%
after = BuCC.d11BCO2().create(numel(co2s));
after.collate("species_calibration").collate("d11B_measured").assignToAll("value",perturbation_d11B4_minimum);

after.collate("boron").collate("pH").assignToAll("pValue",NaN);
after.collate("boron").assignToAll("epsilon",27.2);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

after.collate("carbonate_chemistry").assignToEach("alkalinity",initial.collate("carbonate_chemistry").collate("alkalinity"));

after.calculate();

%
before_pH = initial.collate("carbonate_chemistry").collate("pH").collate("pValue");
after_pH = after.collate("carbonate_chemistry").collate("pH").collate("pValue");

pH_difference("co2") = after_pH-before_pH;

%% Epsilon
saturation_state_maximum = 10.7; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 21.4; % From bulk temperatures
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = 17; % Horita
mg = 28; % Horita

epsilons = [27.2-(0.6*3):0.2:27.2+(0.6*3)];

%
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4_minimum = nanmin(perturbation_data.d11B);

d11B4_change = background_d11B4-perturbation_d11B4_minimum;

% Process with original bulk temperature for comparison
initial = BuCC.d11BCO2().create(numel(epsilons));
initial.collate("species_calibration").collate("d11B_measured").assignToAll("value",background_d11B4);

initial.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial.collate("boron").assignToEach("epsilon",epsilons);
initial.collate("boron").collate("d11B_sw").assignToAll("value",NaN);

cc = initial.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

initial.collate("carbonate_chemistry").collate("atmospheric_co2").assignToAll("partial_pressure",co2_minimum);
initial.collate("carbonate_chemistry").assignToAll("saturation_state",saturation_state_maximum);

initial.calculate();

%
after = BuCC.d11BCO2().create(numel(epsilons));
after.collate("species_calibration").collate("d11B_measured").assignToAll("value",perturbation_d11B4_minimum);

after.collate("boron").collate("pH").assignToAll("pValue",NaN);
after.collate("boron").assignToEach("epsilon",epsilons);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",pressure);
cc.assignToAll("atmospheric_pressure",1);
cc.assignToAll("calcium",ca);
cc.assignToAll("magnesium",mg);

after.collate("carbonate_chemistry").assignToEach("alkalinity",initial.collate("carbonate_chemistry").collate("alkalinity"));

after.calculate();

%
before_pH = initial.collate("carbonate_chemistry").collate("pH").collate("pValue");
after_pH = after.collate("carbonate_chemistry").collate("pH").collate("pValue");

pH_difference("epsilon") = after_pH-before_pH;

%% Saving
output_file = fopen("./../../Data/Minimum_pH_Variation.json","w");
to_encode = [];
to_encode.calcium = cas;
to_encode.pH_difference = pH_difference("calcium");
json = jsonencode(to_encode);
fprintf(output_file,"["+json+",");

to_encode = [];
to_encode.magnesium = mgs;
to_encode.pH_difference = pH_difference("magnesium");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.temperature_change = temperature_changes;
to_encode.pH_difference = pH_difference("temperature_change");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.initial_temperature = initial_temperatures;
to_encode.pH_difference = pH_difference("initial_temperature");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.saturation_state = saturation_states;
to_encode.pH_difference = pH_difference("saturation_state");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.co2 = co2s;
to_encode.pH_difference = pH_difference("co2");
json = jsonencode(to_encode);
fprintf(output_file,json+",");

to_encode = [];
to_encode.epsilon = epsilons;
to_encode.pH_difference = pH_difference("epsilon");
json = jsonencode(to_encode);
fprintf(output_file,json+"]");

fclose(output_file);

