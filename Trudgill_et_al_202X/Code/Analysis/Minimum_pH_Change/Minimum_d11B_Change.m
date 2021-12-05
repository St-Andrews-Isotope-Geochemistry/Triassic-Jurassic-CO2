
%%
clear

data_directory = "./../../../Data/";
number_of_samples = 10000;

boron_data = readtable(data_directory+"/Boron/TJ_d11B.xlsx","Sheet","Delta_Temperature");
boron_data_sorted = sortrows(boron_data,"d11B");
perturbation_boron_data = boron_data_sorted(1:2,:);

%%
perturbation_boron_bin_width = 0.001;
perturbation_boron_bin_edges = 5:perturbation_boron_bin_width:15;
perturbation_boron_bin_midpoints = perturbation_boron_bin_edges(1:end-1)+perturbation_boron_bin_width/2;

for sampler_index = 1:height(perturbation_boron_data)
    perturbation_boron_samplers(sampler_index) = Geochemistry_Helpers.Sampler(perturbation_boron_bin_edges,"Gaussian",[perturbation_boron_data.d11B(sampler_index),perturbation_boron_data.d11B_uncertainty(sampler_index)],"latin_hypercube").normalise();
end

combined_perturbation_boron_distribution = Geochemistry_Helpers.Distribution(perturbation_boron_bin_edges,"Manual",prod(perturbation_boron_samplers.probabilities)).normalise();

%% Minimum Change
maximum_perturbation_boron = combined_perturbation_boron_distribution.quantile(0.95);
maximum_perturbation_boron_uncertainty = combined_perturbation_boron_distribution.standard_deviation();

minimum_pH_change_parameters = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Input.json"));
minimum_initial_boron = minimum_pH_change_parameters.initial_d11B;

minimum_d11B_change = maximum_perturbation_boron(1)-minimum_initial_boron(1);
minimum_d11B_change_uncertainty = maximum_perturbation_boron_uncertainty;

%% Save results
minimum_pH_change_parameters.d11B_change = round([minimum_d11B_change,minimum_d11B_change_uncertainty],3);

writeInputJSON(data_directory+"/Minimum_pH_Change/Input.json",minimum_pH_change_parameters);