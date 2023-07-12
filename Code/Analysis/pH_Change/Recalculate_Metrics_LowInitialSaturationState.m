clear

%% Load data
data_directory = "./../../../Data/";
boron_data = readtable(data_directory+"/Boron/TJ_d11B.xlsx","Sheet","Delta_Temperature");

posterior_file = jsondecode(fileread(data_directory+"/pH_change/posterior.json"));
interpolation_ages = posterior_file.age;
posterior = posterior_file.posterior;

quantile_level = 1:99;
quantiles = (quantile_level/100);

raw_minimum_metrics = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Metrics.json"));

%% Refine by saturation state
valid = [posterior.saturation_state_initial]<5.5;

%%
d11B_measured_evolutions = [posterior(valid).d11B_measured];
pH_evolutions = [posterior(valid).pH];
co2_evolutions = [posterior(valid).co2];
saturation_state_evolutions = [posterior(valid).saturation_state];

delta_saturation_state = [posterior(valid).saturation_state_initial]-min(saturation_state_evolutions);
delta_saturation_state_distribution = Geochemistry_Helpers.Distribution.fromSamples([],delta_saturation_state);

delta_pH = pH_evolutions(end,:)-min(pH_evolutions);
delta_pH_distribution = Geochemistry_Helpers.Distribution.fromSamples([],delta_pH);

%% Get the initial subsample
for distribution_index = 1:size(co2_evolutions(valid),1)
    evolutions.d11B_measured_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],d11B_measured_evolutions(distribution_index,:));
    evolutions.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],pH_evolutions(distribution_index,:));
    evolutions.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],co2_evolutions(distribution_index,:));
    evolutions.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,saturation_state_evolutions(distribution_index,:));
end
evolutions.d11B_measured_quantiles = evolutions.d11B_measured_distributions.quantile(quantiles);
evolutions.pH_quantiles = evolutions.pH_distributions.quantile(quantiles);
evolutions.co2_quantiles = evolutions.co2_distributions.quantile(quantiles);
evolutions.saturation_state_quantiles = evolutions.saturation_state_distributions.quantile(quantiles);


%%
filename = "/pH_Change/TJ_CO2_Evolutions_Metrics_LowInitialSaturationState.csv";

writematrix(evolutions.d11B_measured_quantiles',data_directory+filename,'WriteMode','overwrite');
writematrix(evolutions.pH_quantiles',data_directory+filename,'WriteMode','append');
writematrix(evolutions.saturation_state_quantiles',data_directory+filename,'WriteMode','append');
writematrix(evolutions.co2_quantiles',data_directory+filename,'WriteMode','append');

