clear

%% Load data
data_directory = "./../../../Data/";
boron_data = readtable(data_directory+"/Boron/TJ_d11B.xlsx","Sheet","Delta_Temperature");

minimum_pH_ensemble = readmatrix(data_directory+"/Minimum_pH_Change/TJ_Minimum_pH_Ensemble.csv");

d11B_measured = minimum_pH_ensemble(1:2,:);
pH = minimum_pH_ensemble(3:4,:);

d11B_measured_background = d11B_measured(1,:);
d11B_measured_perturbation = d11B_measured(2,:);
pH_background = pH(1,:);
pH_perturbation = pH(2,:);

% If a value in either array is imaginary then we don't want it
real_pH = imag(pH_perturbation)==0 & imag(pH_background)==0;

d11B_measured_background_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:20,d11B_measured_background(real_pH)).normalise();
d11B_measured_perturbation_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:20,d11B_measured_perturbation(real_pH)).normalise();
pH_background_distribution = Geochemistry_Helpers.Distribution.fromSamples(5:0.01:10,pH_background(real_pH)).normalise();
pH_perturbation_distribution = Geochemistry_Helpers.Distribution.fromSamples(5:0.01:10,pH_perturbation(real_pH)).normalise();

quantile_level = 1:99;
quantiles = (quantile_level/100);




d11B_measured_background_quantiles = d11B_measured_background_distribution.quantile(quantiles);
d11B_measured_perturbation_quantiles = d11B_measured_perturbation_distribution.quantile(quantiles);

pH_background_quantiles = pH_background_distribution.quantile(quantiles);
pH_perturbation_quantiles = pH_perturbation_distribution.quantile(quantiles);


%%
filename = "/Minimum_pH_Change/Quantiles.csv";

writematrix(d11B_measured_background_quantiles,data_directory+filename,'WriteMode','overwrite');
writematrix(d11B_measured_perturbation_quantiles,data_directory+filename,'WriteMode','append');

writematrix(pH_background_quantiles,data_directory+filename,'WriteMode','append');
writematrix(pH_perturbation_quantiles,data_directory+filename,'WriteMode','append');

