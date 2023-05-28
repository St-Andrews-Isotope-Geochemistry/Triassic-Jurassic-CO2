% Generate temperature distributions which are used to do the pH calculations

%%
clear

data_directory = "./../../../Data/";
number_of_samples = 10000;

import_options = detectImportOptions(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx",'NumHeaderLines',2,"Sheet","Temperature_Calibrations");
import_options.VariableNamesRange = "A1:R1";
boron_data = readtable(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx",import_options);
background_data = boron_data(boron_data.age>=boron_data.age(9),:);

%%
temperature_bin_width = 0.1;
temperature_bin_edges = 10:0.1:60;
temperature_bin_midpoints = temperature_bin_edges(1:end-1)+temperature_bin_width/2;

for index = 1:height(background_data)
    if ~isnan(background_data.oneil_temperature(index))
        temperature_samplers(index) = Geochemistry_Helpers.Sampler(temperature_bin_edges,"Gaussian",[background_data.oneil_temperature(index),background_data.oneil_temperature_uncertainty(index)],"latin_hypercube").normalise();
        temperature_samplers(index).getSamples(number_of_samples).shuffle();
        temperature_samples(:,index) = temperature_samplers(index).samples;
    end
end
oneil_temperature_means = mean(temperature_samples,2);
oneil_temperature_mean_distribution = Geochemistry_Helpers.Distribution.fromSamples(temperature_bin_edges,oneil_temperature_means).normalise();

for index = 1:height(background_data)
    if ~isnan(background_data.oneil_temperature(index))
        temperature_samplers(index) = Geochemistry_Helpers.Sampler(temperature_bin_edges,"Gaussian",[background_data.kim_oneil_temperature(index),background_data.kim_oneil_temperature_uncertainty(index)],"latin_hypercube").normalise();
        temperature_samplers(index).getSamples(number_of_samples).shuffle();
        temperature_samples(:,index) = temperature_samplers(index).samples;
    end
end
kim_oneil_temperature_means = mean(temperature_samples,2);
kim_oneil_temperature_mean_distribution = Geochemistry_Helpers.Distribution.fromSamples(temperature_bin_edges,kim_oneil_temperature_means).normalise();

for index = 1:height(background_data)
    if ~isnan(background_data.oneil_temperature(index))
        temperature_samplers(index) = Geochemistry_Helpers.Sampler(temperature_bin_edges,"Gaussian",[background_data.hansen_temperature(index),background_data.hansen_temperature_uncertainty(index)],"latin_hypercube").normalise();
        temperature_samplers(index).getSamples(number_of_samples).shuffle();
        temperature_samples(:,index) = temperature_samplers(index).samples;
    end
end
hansen_temperature_means = mean(temperature_samples,2);
hansen_temperature_mean_distribution = Geochemistry_Helpers.Distribution.fromSamples(temperature_bin_edges,hansen_temperature_means).normalise();

for index = 1:height(background_data)
    if ~isnan(background_data.oneil_temperature(index))
        temperature_samplers(index) = Geochemistry_Helpers.Sampler(temperature_bin_edges,"Gaussian",[background_data.anderson_arthur_temperature(index),background_data.anderson_arthur_temperature_uncertainty(index)],"latin_hypercube").normalise();
        temperature_samplers(index).getSamples(number_of_samples).shuffle();
        temperature_samples(:,index) = temperature_samplers(index).samples;
    end
end
anderson_arthur_temperature_means = mean(temperature_samples,2);
anderson_arthur_temperature_mean_distribution = Geochemistry_Helpers.Distribution.fromSamples(temperature_bin_edges,anderson_arthur_temperature_means).normalise();

combined_probabilities = [oneil_temperature_mean_distribution.probabilities,kim_oneil_temperature_mean_distribution.probabilities,hansen_temperature_mean_distribution.probabilities,anderson_arthur_temperature_mean_distribution.probabilities];
combined_temperature_distribution = Geochemistry_Helpers.Distribution(temperature_bin_edges,"Manual",sum(combined_probabilities,2)).normalise();

%% Fit
temperature_fit_function = @(a1,b1,c1,a2,b2,c2,x) (a1/(c1*sqrt(2*pi)))*exp(-0.5*((x-b1)/c1).^2) + (a2/(c2*sqrt(2*pi)))*exp(-0.5*((x-b2)/c2).^2);
temperature_fit = fit(temperature_bin_midpoints',combined_temperature_distribution.probabilities,temperature_fit_function,'StartPoint',[0.05,18,1,0.05,21,1]);

%% Validate the fit
fit_output = temperature_fit(temperature_bin_edges);

test_distribution(1) = Geochemistry_Helpers.Distribution(temperature_bin_edges,"Gaussian",[temperature_fit.b1,temperature_fit.c1]);
test_distribution(2) = Geochemistry_Helpers.Distribution(temperature_bin_edges,"Gaussian",[temperature_fit.b2,temperature_fit.c2]);

test_distribution(1).probabilities = test_distribution(1).probabilities*temperature_fit.a1;
test_distribution(2).probabilities = test_distribution(2).probabilities*temperature_fit.a2;

test_distribution(3) = Geochemistry_Helpers.Distribution(temperature_bin_edges,"Manual",test_distribution(1).probabilities+test_distribution(2).probabilities);

%
clf
figure(1)
hold on
combined_temperature_distribution.plot();
plot(temperature_bin_edges,fit_output);
test_distribution(3).plot();

%% Save results
pH_change_parameters = jsondecode(fileread(data_directory+"/pH_Change/Input.json"));
pH_change_parameters.initial_temperature = round([temperature_fit.a1,temperature_fit.b1,temperature_fit.c1,temperature_fit.a2,temperature_fit.b2,temperature_fit.c2],3);

writeInputJSON(data_directory+"/pH_Change/Input.json",pH_change_parameters);
