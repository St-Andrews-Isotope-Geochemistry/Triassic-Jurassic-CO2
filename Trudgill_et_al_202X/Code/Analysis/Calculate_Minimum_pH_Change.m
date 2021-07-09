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
data_directory = "./../../Data/";

raw_input_parameters = fileread(data_directory+"Minimum_pH_Change_Input.json");
input_parameters = jsondecode(raw_input_parameters);

%% Set up distributions and sample them
% Choose the number of statistical samples
number_of_samples = 10000;

% Set up a distribution for each variable
% Uncertainties should be tuned to represent 2sigma uncertainty
% If not known, be generous
initial_co2_distribution = Geochemistry_Helpers.Distribution(10e-6:10e-6:1000e-6,"Gaussian",input_parameters.co2/1e6).normalise();
initial_saturation_distribution = Geochemistry_Helpers.Distribution(6:0.01:15,"Gaussian",input_parameters.saturation_state).normalise();
initial_temperature_distribution = Geochemistry_Helpers.Distribution(0:0.01:30,"Gaussian",input_parameters.initial_temperature).normalise();
temperature_change_distribution = Geochemistry_Helpers.Distribution(0:0.01:10,"Gaussian",input_parameters.temperature_change).normalise();
ca_distribution = Geochemistry_Helpers.Distribution(0:0.01:20,"Flat",input_parameters.calcium).normalise();
mg_distribution = Geochemistry_Helpers.Distribution(20:0.01:61,"Flat",input_parameters.magnesium).normalise();
epsilon_distribution = Geochemistry_Helpers.Distribution(27:0.01:31,"Gaussian",input_parameters.epsilon).normalise();
initial_d11B_distribution = Geochemistry_Helpers.Distribution(0:0.01:40,"Gaussian",input_parameters.initial_d11B).normalise();
d11B_change_distribution = Geochemistry_Helpers.Distribution(-8:0.01:0,"Gaussian",input_parameters.d11B_change).normalise();

% Create samplers for those distributions
initial_co2_sampler = Geochemistry_Helpers.Sampler(initial_co2_distribution,"latin_hypercube");
initial_saturation_sampler = Geochemistry_Helpers.Sampler(initial_saturation_distribution,"latin_hypercube");
initial_temperature_sampler = Geochemistry_Helpers.Sampler(initial_temperature_distribution,"latin_hypercube");
temperature_change_sampler = Geochemistry_Helpers.Sampler(temperature_change_distribution,"latin_hypercube");
ca_sampler = Geochemistry_Helpers.Sampler(ca_distribution,"latin_hypercube");
mg_sampler = Geochemistry_Helpers.Sampler(mg_distribution,"latin_hypercube");
epsilon_sampler = Geochemistry_Helpers.Sampler(epsilon_distribution,"latin_hypercube");
initial_d11B_sampler = Geochemistry_Helpers.Sampler(initial_d11B_distribution,"latin_hypercube");
d11B_change_sampler = Geochemistry_Helpers.Sampler(d11B_change_distribution,"latin_hypercube");

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
initial = BuCC.d11BCO2().create(number_of_samples);

% Assign all the necessary values
initial.species_calibration.d11B_measured.assignToEach("value",initial_d11B_sampler.samples);

initial.boron.pH.assignToAll("pValue",NaN);
initial.boron.d11B_sw.assignToAll("value",NaN);
initial.boron.assignToEach("epsilon",epsilon_sampler.samples);

cc = initial.carbonate_chemistry;
cc.assignToAll("units"," mol/kg");
cc.assignToEach("temperature",initial_temperature_sampler.samples);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",oceanic_pressure);
cc.assignToAll("atmospheric_pressure",atmospheric_pressure);
cc.assignToEach("calcium",ca_sampler.samples);
cc.assignToEach("magnesium",mg_sampler.samples);

initial.carbonate_chemistry.atmospheric_co2.assignToEach("partial_pressure",initial_co2_sampler.samples);
initial.carbonate_chemistry.assignToEach("saturation_state",initial_saturation_sampler.samples);

myami = MyAMI.MyAMI("Precalculated",true);
initial.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

% Do carbonate chemistry calculations
initial.calculate();

%% Process minimum d11B measured
% Set up an array to do calculations
after = BuCC.d11BCO2().create(number_of_samples);

% Assign all the necessary variables
after.species_calibration.d11B_measured.assignToEach("value",initial_d11B_sampler.samples+d11B_change_sampler.samples);

after.boron.pH.assignToAll("pValue",NaN);
after.boron.d11B_sw.assignToEach("value",initial.collate("boron").collate("d11B_sw").collate("value"));
after.boron.assignToEach("epsilon",epsilon_sampler.samples);

cc = after.carbonate_chemistry;
cc.assignToAll("units"," mol/kg");
cc.assignToEach("temperature",initial_temperature_sampler.samples+temperature_change_sampler.samples);
cc.assignToAll("salinity",salinity);
cc.assignToAll("oceanic_pressure",oceanic_pressure);
cc.assignToAll("atmospheric_pressure",atmospheric_pressure);
cc.assignToEach("calcium",ca_sampler.samples);
cc.assignToEach("magnesium",mg_sampler.samples);

after.carbonate_chemistry.assignToEach("alkalinity",initial.collate("carbonate_chemistry").collate("alkalinity"));
after.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

% Do carbonate chemistry calculations
after.calculate();

%% Postprocessing
% Collate the information
pH_after = after.carbonate_chemistry.pH.pValue;
pH_initial = initial.carbonate_chemistry.pH.pValue;

% If a value in either array is imaginary then we don't want it
real_pH = imag(pH_after)==0 & imag(pH_initial)==0;

% Create distributions from the samples
pH_after_distribution = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,pH_after(real_pH)).normalise();
pH_initial_distribution = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,pH_initial(real_pH)).normalise();

pH_change = pH_after(real_pH) - pH_initial(real_pH);
pH_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(-1:0.01:1,pH_change).normalise();

%% Statistics
% Can also calculate metrics to describe the distribution
pH_change_99 = pH_change_distribution.quantile(0.99);
pH_change_95 = pH_change_distribution.quantile(0.95);
pH_change_50 = pH_change_distribution.quantile(0.5);

pH_initial_99 = pH_initial_distribution.quantile(0.99);
pH_initial_95 = pH_initial_distribution.quantile(0.95);
pH_initial_50 = pH_initial_distribution.quantile(0.5);

pH_after_99 = pH_after_distribution.quantile(0.99);
pH_after_95 = pH_after_distribution.quantile(0.95);
pH_after_50 = pH_after_distribution.quantile(0.5);

%% Save
data_directory = "./../../Data/";
filename = "TJ_Minimum_pH_Ensemble.csv";

delete(data_directory+filename);

writematrix([initial_d11B_sampler.samples;initial_d11B_sampler.samples+d11B_change_sampler.samples],data_directory+filename,'WriteMode','append');
writematrix([pH_initial,pH_after]',data_directory+filename,'WriteMode','append');
writematrix([initial_temperature_sampler.samples;initial_temperature_sampler.samples+temperature_change_sampler.samples],data_directory+filename,'WriteMode','append');

writematrix(initial_co2_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(initial_saturation_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(ca_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(mg_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(epsilon_sampler.samples,data_directory+filename,'WriteMode','append');
writematrix(initial.boron.d11B_sw.value',data_directory+filename,'WriteMode','append');

filename = "Minimum_pH_Change_Metrics.json";

fileID = fopen(data_directory+filename,"w");
fwrite(fileID,"{"+newline+string(char(9))+'"initial":')
fwrite(fileID,"{"+newline+string(char(9))+string(char(9))+'"pH_median":'+pH_initial_50+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_95":'+pH_initial_95+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_99":'+pH_initial_99+newline+string(char(9))+"},"+newline);
fwrite(fileID,string(char(9))+'"after":')
fwrite(fileID,"{"+newline+string(char(9))+string(char(9))+'"pH_median":'+pH_after_50+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_95":'+pH_after_95+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_99":'+pH_after_99+newline+string(char(9))+"},"+newline);
fwrite(fileID,string(char(9))+'"change":')
fwrite(fileID,"{"+newline+string(char(9))+string(char(9))+'"pH_median":'+pH_change_50+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_95":'+pH_change_95+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_99":'+pH_change_99+newline+string(char(9))+"}"+newline);
fwrite(fileID,"}");
fclose(fileID);