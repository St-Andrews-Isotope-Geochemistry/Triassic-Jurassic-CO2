
%%
clear

data_directory = "./../../../Data/";
number_of_samples = 10000;

boron_data = readtable(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx","Sheet","Temperature_Calibrations");
boron_data = boron_data((~boron_data.diagenetic_alteration & ~boron_data.al_ca_reject),:);

boron_data_sorted = sortrows(boron_data,"d11B");

delta_temperature_bin_edges = -20:0.1:20;

%%
perturbation.samples = boron_data(boron_data.age>=min(boron_data_sorted.age(1:2)) & boron_data.age<=max(boron_data_sorted.age(1:2)),:);

for sampler_index = 1:height(perturbation.samples)
    perturbation.samplers(sampler_index) = Geochemistry_Helpers.Sampler(delta_temperature_bin_edges,"Gaussian",[perturbation.samples.delta_temperature(sampler_index),perturbation.samples.hansen_temperature_uncertainty(sampler_index)],'latin_hypercube').normalise();
end

perturbation.combined_sampler = Geochemistry_Helpers.Sampler(delta_temperature_bin_edges',"Manual",prod(perturbation.samplers.probabilities),"latin_hyercube").normalise();

minimum_temperature_rise = perturbation.combined_sampler.quantile(0.05);
minimum_temperature_rise_uncertainty = perturbation.combined_sampler.standard_deviation();

average_temperature_rise = perturbation.combined_sampler.quantile(0.5);
average_temperature_rise_uncertainty = perturbation.combined_sampler.standard_deviation()*2;

%% Save results
minimum_pH_change_parameters = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Input.json"));
minimum_pH_change_parameters.temperature_change = round([minimum_temperature_rise,minimum_temperature_rise_uncertainty],3);

writeInputJSON(data_directory+"/Minimum_pH_Change/Input.json",minimum_pH_change_parameters);