clear

%% Load data
data_directory = "./../../../Data/";
boron_data = readtable(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx","Sheet","Temperature_Calibrations");

% prior_file = jsondecode(fileread(data_directory+"/pH_change/prior_standard_smooth.json"));
% intermediate_file = jsondecode(fileread(data_directory+"/pH_change/intermediate_standard_smooth.json"));
posterior_file = jsondecode(fileread(data_directory+"/pH_change/posterior.json"));
interpolation_ages = posterior_file.age;

posterior = posterior_file.posterior;

%% Extract relevant parameters for posterior
d11B_measured_evolutions = [posterior.d11B_measured];
pH_evolutions = [posterior.pH];
co2_evolutions = [posterior.co2];
saturation_state_evolutions = [posterior.saturation_state];

delta_pH = pH_evolutions-pH_evolutions(end,:);
delta_co2 = log2(co2_evolutions)-log2(co2_evolutions(end,:));

%% Get the summary stats
quantiles = [0.025,0.05,0.5,0.095,0.975];
for distribution_index = 1:size(co2_evolutions,1)
    evolutions.d11B_measured_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],d11B_measured_evolutions(distribution_index,:));
    evolutions.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],pH_evolutions(distribution_index,:));
    evolutions.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],co2_evolutions(distribution_index,:));
    evolutions.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],saturation_state_evolutions(distribution_index,:));
    evolutions.delta_co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],delta_co2(distribution_index,:));
end
evolutions.d11B_measured_quantiles = evolutions.d11B_measured_distributions.quantile(quantiles);
evolutions.pH_quantiles = evolutions.pH_distributions.quantile(quantiles);
evolutions.co2_quantiles = evolutions.co2_distributions.quantile(quantiles);
evolutions.saturation_state_quantiles = evolutions.saturation_state_distributions.quantile(quantiles);
evolutions.delta_co2_quantiles = evolutions.delta_co2_distributions.quantile(quantiles);

%% Make tables
d11B_table = table();
d11B_table.age = interpolation_ages;
d11B_table.quantile_0025 = evolutions.d11B_measured_quantiles(:,1);
d11B_table.quantile_0050 = evolutions.d11B_measured_quantiles(:,2);
d11B_table.quantile_0500 = evolutions.d11B_measured_quantiles(:,3);
d11B_table.quantile_0950 = evolutions.d11B_measured_quantiles(:,4);
d11B_table.quantile_0975 = evolutions.d11B_measured_quantiles(:,5);

pH_table = table();
pH_table.age = interpolation_ages;
pH_table.quantile_0025 = evolutions.pH_quantiles(:,1);
pH_table.quantile_0050 = evolutions.pH_quantiles(:,2);
pH_table.quantile_0500 = evolutions.pH_quantiles(:,3);
pH_table.quantile_0950 = evolutions.pH_quantiles(:,4);
pH_table.quantile_0975 = evolutions.pH_quantiles(:,5);

co2_table = table();
co2_table.age = interpolation_ages;
co2_table.quantile_0025 = evolutions.co2_quantiles(:,1);
co2_table.quantile_0050 = evolutions.co2_quantiles(:,2);
co2_table.quantile_0500 = evolutions.co2_quantiles(:,3);
co2_table.quantile_0950 = evolutions.co2_quantiles(:,4);
co2_table.quantile_0975 = evolutions.co2_quantiles(:,5);

saturation_state_table = table();
saturation_state_table.age = interpolation_ages;
saturation_state_table.quantile_0025 = evolutions.saturation_state_quantiles(:,1);
saturation_state_table.quantile_0050 = evolutions.saturation_state_quantiles(:,2);
saturation_state_table.quantile_0500 = evolutions.saturation_state_quantiles(:,3);
saturation_state_table.quantile_0950 = evolutions.saturation_state_quantiles(:,4);
saturation_state_table.quantile_0975 = evolutions.saturation_state_quantiles(:,5);

delta_co2_table = table();
delta_co2_table.age = interpolation_ages;
delta_co2_table.quantile_0025 = evolutions.delta_co2_quantiles(:,1);
delta_co2_table.quantile_0050 = evolutions.delta_co2_quantiles(:,2);
delta_co2_table.quantile_0500 = evolutions.delta_co2_quantiles(:,3);
delta_co2_table.quantile_0950 = evolutions.delta_co2_quantiles(:,4);
delta_co2_table.quantile_0975 = evolutions.delta_co2_quantiles(:,5);

%% Save to file
filename = data_directory+"/pH_Change/posterior_standard_smooth_metrics.xlsx";

writetable(d11B_table,filename,'Sheet','d11B_measured');
writetable(pH_table,filename,'Sheet','pH');
writetable(co2_table,filename,'Sheet','CO2');
writetable(delta_co2_table,filename,'Sheet','Delta_CO2');
writetable(saturation_state_table,filename,'Sheet','Saturation_State');