%% Minimum pH Change


%%
clear

%%
temperature_data = readtable("./../../Data/TJ_Temperature.xlsx");
boron_data = readtable("./../../Data/TJ_d11B_pH.xlsx");

saturation_state_maximum = 10.7; % From Ridgwell modelling
co2_minimum = 500e-6; % From Witkowski

temperature = 21.4; % From bulk temperatures
salinity = 35; % Assumed
pressure = 0; % Assumed?
ca = 17; % Horita
mg = 28; % Horita

%% Need to interpolate temperature onto the samples
boron_data.temperature = interp1(temperature_data.depth,temperature_data.temperature,boron_data.height);
boron_data.temperature_uncertainty = interp1(temperature_data.depth,temperature_data.temperature_uncertainty,boron_data.height);

%%
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);
background_temperature = nanmean(background_data.temperature);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4_minimum = nanmin(perturbation_data.d11B);
perturbation_temperature = nanmin(perturbation_data.temperature);

d11B4_change = background_d11B4-perturbation_d11B4_minimum;

%% Process with original bulk temperature for comparison
initial = BuCC.d11BCO2();
initial.species_calibration.d11B_measured.value = background_d11B4;

initial.boron.pH.pValue = NaN;
initial.boron.epsilon = 27.2+2*0.06;
initial.boron.d11B_sw.value = NaN;

background_conditions = BuCC.Conditions("temperature",temperature,"salinity",salinity,"oceanic_pressure",pressure,"calcium",ca,"magnesium",mg,"atmospheric_pressure",1);
initial.carbonate_chemistry.setConditions(background_conditions);

initial.carbonate_chemistry.atmospheric_co2.partial_pressure = co2_minimum;
initial.carbonate_chemistry.saturation_state = saturation_state_maximum;

initial.calculate();

%% Now imagine changing d11B_4 (no temperature change)
original_d11B4_change = BuCC.d11BCO2();
original_d11B4_change.species_calibration.d11B_measured.value = perturbation_d11B4_minimum;

original_d11B4_change.boron.pH.pValue = NaN;
original_d11B4_change.boron.epsilon = 27.2+2*0.06;
original_d11B4_change.boron.d11B_sw.value = initial.boron.d11B_sw.value;

background_conditions = BuCC.Conditions("temperature",temperature,"salinity",salinity,"oceanic_pressure",pressure,"calcium",ca,"magnesium",mg,"atmospheric_pressure",1);
original_d11B4_change.carbonate_chemistry.setConditions(background_conditions);

original_d11B4_change.carbonate_chemistry.alkalinity = initial.carbonate_chemistry.alkalinity;

original_d11B4_change.calculate();

%% Instead process with the d18O derived temperature
background = BuCC.d11BCO2();
background.species_calibration.d11B_measured.value = background_d11B4;

background.boron.pH.pValue = NaN;
background.boron.epsilon = 27.2+2*0.06;
background.boron.d11B_sw.value = NaN;

background_conditions = BuCC.Conditions("temperature",background_temperature,"salinity",salinity,"oceanic_pressure",pressure,"calcium",ca,"magnesium",mg,"atmospheric_pressure",1);
background.carbonate_chemistry.setConditions(background_conditions);

background.carbonate_chemistry.atmospheric_co2.partial_pressure = co2_minimum;
background.carbonate_chemistry.saturation_state = saturation_state_maximum;

background.calculate();

%% Now imagine changing d11B_4 (no temperature change)
perturbation_d11B4_change = BuCC.d11BCO2();
perturbation_d11B4_change.species_calibration.d11B_measured.value = perturbation_d11B4_minimum;

perturbation_d11B4_change.boron.pH.pValue = NaN;
perturbation_d11B4_change.boron.epsilon = 27.2+2*0.06;
perturbation_d11B4_change.boron.d11B_sw.value = background.boron.d11B_sw.value;

perturbation_d11B4_change.carbonate_chemistry.setConditions(background_conditions);
perturbation_d11B4_change.carbonate_chemistry.alkalinity = background.carbonate_chemistry.alkalinity;

perturbation_d11B4_change.calculate();

%% Now change both d11B_4 and temperature
perturbation_d11B4_temperature_change = BuCC.d11BCO2();
perturbation_d11B4_temperature_change.species_calibration.d11B_measured.value = perturbation_d11B4_minimum;

perturbation_d11B4_temperature_change.boron.pH.pValue = NaN;
perturbation_d11B4_temperature_change.boron.epsilon = 27.2+2*0.06;
perturbation_d11B4_temperature_change.boron.d11B_sw.value = background.boron.d11B_sw.value;

perturbation_conditions = BuCC.Conditions("temperature",perturbation_temperature,"salinity",salinity,"oceanic_pressure",pressure,"calcium",ca,"magnesium",mg,"atmospheric_pressure",1);
perturbation_d11B4_temperature_change.carbonate_chemistry.setConditions(perturbation_conditions);

perturbation_d11B4_temperature_change.carbonate_chemistry.alkalinity = background.carbonate_chemistry.alkalinity;

perturbation_d11B4_temperature_change.calculate();

%% Collate objects and values
objects = [initial,original_d11B4_change,background,perturbation_d11B4_change,perturbation_d11B4_temperature_change];

pHs = objects.collate("carbonate_chemistry").collate("pH").collate("pValue");
co2s = objects.collate("carbonate_chemistry").collate("atmospheric_co2").collate("partial_pressure")*1e6;
alkalinities = objects.collate("carbonate_chemistry").collate("alkalinity")*1e6;
