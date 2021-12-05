
%%
clear

data_directory = "./../../../Data/";
number_of_samples = 10000;

boron_data = readtable(data_directory+"/Boron/TJ_d11B.xlsx","Sheet","Delta_Temperature");
background_boron_data = boron_data(1:9,:);

%%
background_boron_bin_width = 0.001;
background_boron_bin_edges = 5:background_boron_bin_width:15;
background_boron_bin_midpoints = background_boron_bin_edges(1:end-1)+background_boron_bin_width/2;

for sampler_index = 1:height(background_boron_data)
    background_boron_samplers(sampler_index) = Geochemistry_Helpers.Sampler(background_boron_bin_edges,"Gaussian",[background_boron_data.d11B(sampler_index),background_boron_data.d11B_uncertainty(sampler_index)],"latin_hypercube").normalise();
end

combined_initial_boron_distribution = Geochemistry_Helpers.Distribution(background_boron_bin_edges,"Manual",prod(reshape(background_boron_samplers.probabilities,[],numel(background_boron_bin_midpoints)))).normalise();

%% Maximum temperature
minimum_boron = combined_initial_boron_distribution.quantile(0.05);
minimum_boron_uncertainty = combined_initial_boron_distribution.standard_deviation();

%% Save results
minimum_pH_change_parameters = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Input.json"));
minimum_pH_change_parameters.initial_d11B = round([minimum_boron,minimum_boron_uncertainty],3);

writeInputJSON(data_directory+"/Minimum_pH_Change/Input.json",minimum_pH_change_parameters);