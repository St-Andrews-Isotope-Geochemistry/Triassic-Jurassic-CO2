% Minimum pH change sampler
% There are uncertainties on a lot of parameters necessary to calculate the
% change in pH
% To account for those uncertainties we want distributions for static parameters:
% - initial temperature
% - initial saturation state
% - initial CO2
% - ca
% - mg
% - epsilon
%
% And distributions for change parameters:
% - d11B change
% - temperature change

% We'll work at 95% uncertainty (2sigma)
clear

%% Load in data
data_directory = "./../../../Data/";

input_parameters = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Input.json"));

%% Set up distributions and sample them
% Choose the number of statistical samples
number_of_samples = 10000;

% Set up a sampler for each variable
% Uncertainties should be tuned to represent 2sigma uncertainty
% If not known, be generous
initial_co2_sampler = Geochemistry_Helpers.Sampler(10e-6:10e-6:1000e-6,"Gaussian",input_parameters.co2/1e6,"latin_hypercube").normalise();
initial_saturation_sampler = Geochemistry_Helpers.Sampler(6:0.01:15,"Gaussian",input_parameters.saturation_state,"latin_hypercube").normalise();
initial_temperature_sampler = Geochemistry_Helpers.Sampler(0:0.01:30,"Gaussian",input_parameters.initial_temperature,"latin_hypercube").normalise();
temperature_change_sampler = Geochemistry_Helpers.Sampler(0:0.01:10,"Gaussian",input_parameters.temperature_change,"latin_hypercube").normalise();
ca_sampler = Geochemistry_Helpers.Sampler(0:0.01:20,"Flat",input_parameters.calcium,"latin_hypercube").normalise();
mg_sampler = Geochemistry_Helpers.Sampler(20:0.01:61,"Flat",input_parameters.magnesium,"latin_hypercube").normalise();
epsilon_sampler = Geochemistry_Helpers.Sampler(27:0.01:31,"Gaussian",input_parameters.epsilon,"latin_hypercube").normalise();
initial_d11B_sampler = Geochemistry_Helpers.Sampler(0:0.01:40,"Gaussian",input_parameters.initial_d11B,"latin_hypercube").normalise();
d11B_change_sampler = Geochemistry_Helpers.Sampler(-8:0.01:0,"Gaussian",input_parameters.d11B_change,"latin_hypercube").normalise();

% Get samples
initial_co2_sampler.getSamples(number_of_samples).shuffle();
initial_saturation_sampler.getSamples(number_of_samples).shuffle();
initial_temperature_sampler.getSamples(number_of_samples).shuffle();
temperature_change_sampler.getSamples(number_of_samples).shuffle();
ca_sampler.getSamples(number_of_samples).shuffle();
mg_sampler.getSamples(number_of_samples).shuffle();
epsilon_sampler.getSamples(number_of_samples).shuffle();
initial_d11B_sampler.getSamples(number_of_samples).shuffle();
d11B_change_sampler.getSamples(number_of_samples).shuffle();

% Additional parameters
oceanic_pressure = input_parameters.oceanic_pressure;
atmospheric_pressure = input_parameters.atmospheric_pressure;
salinity = input_parameters.salinity;

%% Process background
% Set up an array to do the d11B->CO2 calculations
background.samples = BuCC.d11BCO2().create(number_of_samples);

% Assign all the necessary values
% Species calibration
background.samples.species_calibration.d11B_measured.assignToEach("value",initial_d11B_sampler.samples);

% Boron
background.samples.boron.pH.assignToAll("pValue",NaN);
background.samples.boron.d11B_sw.assignToAll("value",NaN);
background.samples.boron.assignToEach("epsilon",epsilon_sampler.samples);

% Carbonate Chemistry
background.samples.carbonate_chemistry.assignToAll("units"," mol/kg");
background.samples.carbonate_chemistry.assignToEach("temperature",initial_temperature_sampler.samples);
background.samples.carbonate_chemistry.assignToAll("salinity",salinity);
background.samples.carbonate_chemistry.assignToAll("oceanic_pressure",oceanic_pressure);
background.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",atmospheric_pressure);
background.samples.carbonate_chemistry.assignToEach("calcium",ca_sampler.samples);
background.samples.carbonate_chemistry.assignToEach("magnesium",mg_sampler.samples);

background.samples.carbonate_chemistry.atmospheric_co2.assignToEach("partial_pressure",initial_co2_sampler.samples);
background.samples.carbonate_chemistry.assignToEach("saturation_state",initial_saturation_sampler.samples);

myami = MyAMI.MyAMI("Precalculated",true);
background.samples.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

% Do calculations
background.samples.calculate();

%% Process perturbed d11B (minimum)
% Set up an array to do calculations
perturbation.samples = BuCC.d11BCO2().create(number_of_samples);

% Assign all the necessary variables
% Species calibration
perturbation.samples.species_calibration.d11B_measured.assignToEach("value",initial_d11B_sampler.samples+d11B_change_sampler.samples);

% Boron
perturbation.samples.boron.pH.assignToAll("pValue",NaN);
perturbation.samples.boron.d11B_sw.assignToEach("value",background.samples.collate("boron").collate("d11B_sw").collate("value"));
perturbation.samples.boron.assignToEach("epsilon",epsilon_sampler.samples);

% Carbonate Chemistry
perturbation.samples.carbonate_chemistry.assignToAll("units"," mol/kg");
perturbation.samples.carbonate_chemistry.assignToEach("temperature",initial_temperature_sampler.samples+temperature_change_sampler.samples);
perturbation.samples.carbonate_chemistry.assignToAll("salinity",salinity);
perturbation.samples.carbonate_chemistry.assignToAll("oceanic_pressure",oceanic_pressure);
perturbation.samples.carbonate_chemistry.assignToAll("atmospheric_pressure",atmospheric_pressure);
perturbation.samples.carbonate_chemistry.assignToEach("calcium",ca_sampler.samples);
perturbation.samples.carbonate_chemistry.assignToEach("magnesium",mg_sampler.samples);

perturbation.samples.carbonate_chemistry.assignToEach("alkalinity",background.samples.collate("carbonate_chemistry").collate("alkalinity"));
perturbation.samples.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

% Do calculations
perturbation.samples.calculate();

%% Postprocessing
% Collate the information
perturbation.pH = perturbation.samples.carbonate_chemistry.pH.pValue;
background.pH = background.samples.carbonate_chemistry.pH.pValue;

% If a value in either array is imaginary then we don't want it (means we've pushed the value past what is possible on the d11B curve)
real_pH = imag(perturbation.pH)==0 & imag(background.pH )==0;
change.pH = perturbation.pH(real_pH) - background.pH(real_pH);

% Create distributions from the samples
background.pH_distribution = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,background.pH(real_pH)).normalise();
perturbation.pH_distribution = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,perturbation.pH(real_pH)).normalise();
change.pH_distribution = Geochemistry_Helpers.Distribution.fromSamples(-1:0.01:1,change.pH).normalise();

%% Statistics
% Can also calculate metrics to describe the distribution
background.pH_metrics = background.pH_distribution.quantile([0.5,0.95,0.99]);
perturbation.pH_metrics = perturbation.pH_distribution.quantile([0.5,0.95,0.99]);
change.pH_metrics = change.pH_distribution.quantile([0.5,0.95,0.99]);

%% Save
% Full ensemble
data_directory = "./../../../Data/";
filename = "/Minimum_pH_Change/TJ_Minimum_pH_Ensemble.csv";

delete(data_directory+filename);

writematrix([initial_d11B_sampler.samples;initial_d11B_sampler.samples+d11B_change_sampler.samples],data_directory+filename,'WriteMode','append');
writematrix([background.pH,perturbation.pH]',data_directory+filename,'WriteMode','append');
writematrix([initial_temperature_sampler.samples;initial_temperature_sampler.samples+temperature_change_sampler.samples],data_directory+filename,'WriteMode','append');

writematrix(initial_co2_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(initial_saturation_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(ca_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(mg_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(epsilon_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(background.samples.boron.d11B_sw.value',data_directory+filename,'WriteMode','append');

% Metrics
filename = "/Minimum_pH_Change/Metrics.json";

fileID = fopen(data_directory+filename,"w");
fwrite(fileID,"{"+newline+string(char(9))+'"initial":');
fwrite(fileID,"{"+newline+string(char(9))+string(char(9))+'"pH_median":'+background.pH_metrics(1)+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_95":'+background.pH_metrics(2)+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_99":'+background.pH_metrics(3)+newline+string(char(9))+"},"+newline);
fwrite(fileID,string(char(9))+'"after":');
fwrite(fileID,"{"+newline+string(char(9))+string(char(9))+'"pH_median":'+perturbation.pH_metrics(1)+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_95":'+perturbation.pH_metrics(2)+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_99":'+perturbation.pH_metrics(3)+newline+string(char(9))+"},"+newline);
fwrite(fileID,string(char(9))+'"change":');
fwrite(fileID,"{"+newline+string(char(9))+string(char(9))+'"pH_median":'+change.pH_metrics(1)+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_95":'+change.pH_metrics(2)+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_99":'+change.pH_metrics(3)+newline+string(char(9))+"}"+newline);
fwrite(fileID,"}");
fclose(fileID);