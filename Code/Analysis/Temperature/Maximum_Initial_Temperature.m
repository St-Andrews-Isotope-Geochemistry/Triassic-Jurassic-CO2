% Generate temperature distributions which are used to do the pH calculations

%%
clear

data_directory = "./../../../Data/";
number_of_samples = 10000;

boron_data = readtable(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx","Sheet","Temperature_Calibrations");
background_temperature_data = boron_data(boron_data.age>=boron_data.age(9),:);

%%
background_temperature_bin_width = 0.1;
background_temperature_bin_edges = 0:0.1:40;
background_temperature_bin_midpoints = background_temperature_bin_edges(1:end-1)+background_temperature_bin_width/2;

for sampler_index = 1:height(background_temperature_data)
    background_temperature_samplers(sampler_index,1) = Geochemistry_Helpers.Sampler(background_temperature_bin_edges,"Gaussian",[boron_data.oneil_temperature(sampler_index),boron_data.oneil_temperature_uncertainty(sampler_index)],"latin_hypercube").normalise();
    background_temperature_samplers(sampler_index,2) = Geochemistry_Helpers.Sampler(background_temperature_bin_edges,"Gaussian",[boron_data.kim_oneil_temperature(sampler_index),boron_data.kim_oneil_temperature_uncertainty(sampler_index)],"latin_hypercube").normalise();
    background_temperature_samplers(sampler_index,3) = Geochemistry_Helpers.Sampler(background_temperature_bin_edges,"Gaussian",[boron_data.hansen_temperature(sampler_index),boron_data.hansen_temperature_uncertainty(sampler_index)],"latin_hypercube").normalise();
    background_temperature_samplers(sampler_index,4) = Geochemistry_Helpers.Sampler(background_temperature_bin_edges,"Gaussian",[boron_data.anderson_arthur_temperature(sampler_index),boron_data.anderson_arthur_temperature_uncertainty(sampler_index)],"latin_hypercube").normalise();
end

combined_initial_temperature_distribution = Geochemistry_Helpers.Distribution(background_temperature_bin_edges,"Manual",sum(reshape(background_temperature_samplers.probabilities,[],numel(background_temperature_bin_midpoints)),1,"omitnan")).normalise();

%% Maximum temperature
maximum_temperature = combined_initial_temperature_distribution.quantile(0.95);
maximum_temperature_uncertainty = combined_initial_temperature_distribution.standard_deviation();

median_temperature = combined_initial_temperature_distribution.quantile(0.5);
median_temperature_uncertainty = combined_initial_temperature_distribution.standard_deviation()*2;

%% Save results
minimum_pH_change_parameters = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Input.json"));
minimum_pH_change_parameters.initial_temperature = round([maximum_temperature,maximum_temperature_uncertainty],3);

writeInputJSON(data_directory+"/Minimum_pH_Change/Input.json",minimum_pH_change_parameters);