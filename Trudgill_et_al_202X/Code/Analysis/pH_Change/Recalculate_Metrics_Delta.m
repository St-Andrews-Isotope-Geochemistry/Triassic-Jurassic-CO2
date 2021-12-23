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

%%
pH_evolutions = [posterior.pH];
co2_evolutions = [posterior.co2];
saturation_state_evolutions = [posterior.saturation_state];

delta_pH = pH_evolutions-pH_evolutions(end,:);
delta_co2 = log2(co2_evolutions)-log2(co2_evolutions(end,:));
delta_saturation_state = saturation_state_evolutions-saturation_state_evolutions(end,:);

%% Get the initial subsample
for distribution_index = 1:size(delta_pH,1)
    evolutions.delta_pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],delta_pH(distribution_index,:));
    evolutions.delta_co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],delta_co2(distribution_index,:));
    evolutions.delta_saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],delta_saturation_state(distribution_index,:));
end
evolutions.delta_pH_quantiles = evolutions.delta_pH_distributions.quantile(quantiles);
evolutions.delta_co2_quantiles = evolutions.delta_co2_distributions.quantile(quantiles);
evolutions.delta_saturation_state_quantiles = evolutions.delta_saturation_state_distributions.quantile(quantiles);

%%
filename = "/pH_Change/Metrics_Delta.csv";

writematrix(evolutions.delta_pH_quantiles',data_directory+filename,'WriteMode','overwrite');
writematrix(evolutions.delta_saturation_state_quantiles',data_directory+filename,'WriteMode','append');
writematrix(evolutions.delta_co2_quantiles',data_directory+filename,'WriteMode','append');


