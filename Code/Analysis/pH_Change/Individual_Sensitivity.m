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

clear all

%% Load Data
data_directory = "./../../../Data/";
temperature_data = readtable(data_directory+"/Temperature/TJ_d18O_d13C.xlsx");
boron_data = readtable(data_directory+"/Boron/TJ_d11B.xlsx");

% Create a map to hold the results
pH_difference = containers.Map();

%% Ca
saturation_state_maximum = (3.85+11.85)/2; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 20.98; % From bulk temperatures
temperature_change = 4.24;
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = (10+17)/2; % Horita
mg = (28+52)/2; % Horita
epsilon = 27.2;

% For minimum 8-10, 8-17
cas = 13.5-6.5:1:13.5+6.5;

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
initial.collate("boron").assignToAll("epsilon",epsilon);
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
after.collate("boron").assignToAll("epsilon",epsilon);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("temperature",temperature+temperature_change);
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

pH_change = after_pH-before_pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("calcium") = pH_change;

%% Magnesium
saturation_state_maximum = (3.85+11.85)/2; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 20.98; % From bulk temperatures
temperature_change = 4.24;
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = (10+17)/2; % Horita
mg = (28+52)/2; % Horita
epsilon = 27.2;

% For minimum: 50-52, 28-52
mgs = 40-22:1:40+22;

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
initial.collate("boron").assignToAll("epsilon",epsilon);
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
after.collate("boron").assignToAll("epsilon",epsilon);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("temperature",temperature+temperature_change);
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

pH_change = after_pH-before_pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("magnesium") = pH_change;

%% Temperature change
saturation_state_maximum = (3.85+11.85)/2; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 20.98; % From bulk temperatures
temperature_change = 4.24;
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = (10+17)/2; % Horita
mg = (28+52)/2; % Horita
epsilon = 27.2;

% For minimum 2+-0.55 at 1sd, 4.2+-1.1 at 1sd
initial_temperature = temperature;
temperature_changes = 4.24-1.03*1.8:0.1:4.24+1.03*1.8;

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
initial.collate("boron").assignToAll("epsilon",epsilon);
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
after.collate("boron").assignToAll("epsilon",epsilon);
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

pH_change = after_pH-before_pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("temperature_change") = pH_change;

%% Temperature start
saturation_state_maximum = (3.85+11.85)/2; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 20.98; % From bulk temperatures
temperature_change = 4.24;
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = (10+17)/2; % Horita
mg = (28+52)/2; % Horita
epsilon = 27.2;

% For minimum 21.4+-1.65 at 1sd, 14.8+-3.3 at 1sd
initial_temperatures = 20.98-25:1:20.98+25;

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
initial.collate("boron").assignToAll("epsilon",epsilon);
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
after.collate("boron").assignToAll("epsilon",epsilon);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToEach("temperature",initial_temperatures);
cc.assignToAll("temperature",temperature+temperature_change);
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

pH_change = after_pH-before_pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("initial_temperature") = pH_change;

%% Saturation State
saturation_state_maximum = (3.85+11.85)/2; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 20.98; % From bulk temperatures
temperature_change = 4.24;
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = (10+17)/2; % Horita
mg = (28+52)/2; % Horita
epsilon = 27.2;

% For minimum 10.7+-0.15 at 1sd, 5-10.7 for central
saturation_states = 7.85-5.3:0.1:7.85+5.3;

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
initial.collate("boron").assignToAll("epsilon",epsilon);
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
after.collate("boron").assignToAll("epsilon",epsilon);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("temperature",temperature+temperature_change);
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

pH_change = after_pH-before_pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("saturation_state") = pH_change;


%% CO2
saturation_state_maximum = (3.85+11.85)/2; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 20.98; % From bulk temperatures
temperature_change = 4.24;
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = (10+17)/2; % Horita
mg = (28+52)/2; % Horita
epsilon = 27.2;

% For minimum 500e-6 +- 50e-6 at 1sd, 400e-6 - 5000e-6 for central
co2s = 300e-6:50e-6:700e-6;

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
initial.collate("boron").assignToAll("epsilon",epsilon);
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
after.collate("boron").assignToAll("epsilon",epsilon);
after.collate("boron").collate("d11B_sw").assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));

cc = after.collate("carbonate_chemistry");
cc.assignToAll("temperature",temperature);
cc.assignToAll("temperature",temperature+temperature_change);
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

pH_change = after_pH-before_pH;
pH_change(imag(pH_change)~=0) = NaN;

pH_difference("co2") = pH_change;

%% Epsilon
saturation_state_maximum = (3.85+11.85)/2; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 20.98; % From bulk temperatures
temperature_change = 4.24;
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = (10+17)/2; % Horita
mg = (28+52)/2; % Horita
epsilon = 27.2;

% For minimum (27.2+0.3*2) +- 0.15 at 1sd and 27.2 +- 0.3 at 1sd for central
epsilons = 27.2-1.1:0.1:27.2+1.1;

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
cc.assignToAll("temperature",temperature+temperature_change);
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

pH_change = after_pH-before_pH;
pH_change(imag(pH_change)~=0) = NaN;
pH_difference("epsilon") = pH_change;

%% Saving
output_file = fopen(data_directory+"/pH_Change/Individual_Sensitivity.json","w");
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

