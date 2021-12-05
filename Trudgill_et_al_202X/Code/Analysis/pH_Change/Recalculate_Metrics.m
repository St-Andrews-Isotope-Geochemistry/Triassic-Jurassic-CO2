clear all

%% Load data
boron_data = readtable("./../../Data/TJ_d11B.xlsx","Sheet","Delta_Temperature");
evolutions = getValidSamples(getShapedEvolutions("./../../Data/TJ_CO2_Evolutions.csv"));

interpolation_ages = unique(sort([boron_data.age',linspace(min(boron_data.age),max(boron_data.age),80)]));
evolutions.age = repmat(interpolation_ages',1,size(evolutions.pH,2));

quantile_level = 5; % As percent
quantiles = [(quantile_level/100)/2,0.5,1-((quantile_level/100)/2)];
above_quantiles = [0.05,0.5,0.95];

raw_minimum_metrics = jsondecode(fileread("./../../Data/Minimum_pH_Change_Metrics.json"));

%% Get the initial subsample
for distribution_index = 1:size(evolutions.pH,1)
    evolutions.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],evolutions.pH(distribution_index,:));
    evolutions.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],evolutions.co2(distribution_index,:));
    evolutions.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,evolutions.saturation_state(distribution_index,:));
    evolutions.d11B_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:50,evolutions.d11B(distribution_index,:));
    evolutions.alkalinity_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:20000,evolutions.alkalinity(distribution_index,:)*1e6);
end
evolutions.pH_quantiles = evolutions.pH_distributions.quantile(quantiles);
evolutions.saturation_state_quantiles = evolutions.saturation_state_distributions.quantile(quantiles);
evolutions.d11B_quantiles = evolutions.d11B_distributions.quantile(quantiles);
evolutions.co2_quantiles = evolutions.co2_distributions.quantile(quantiles);

evolutions.pH_above_quantiles = evolutions.pH_distributions.quantile(above_quantiles);
evolutions.d11B_above_quantiles = evolutions.d11B_distributions.quantile(above_quantiles);
evolutions.co2_above_quantiles = evolutions.co2_distributions.quantile(above_quantiles);

%% Get the high pH subsample
high_initial_pH_evolutions.boolean = repmat(evolutions.pH(end,:)>round(raw_minimum_metrics.initial.pH_median-0.05,1),size(evolutions.pH,1),1);
high_initial_pH_evolutions.pH = reshape(evolutions.pH(high_initial_pH_evolutions.boolean),size(evolutions.pH,1),[]);
high_initial_pH_evolutions.co2 = reshape(evolutions.co2(high_initial_pH_evolutions.boolean),size(evolutions.pH,1),[]);
high_initial_pH_evolutions.saturation_state = reshape(evolutions.saturation_state(high_initial_pH_evolutions.boolean),size(evolutions.pH,1),[]);
high_initial_pH_evolutions.d11B = reshape(evolutions.d11B(high_initial_pH_evolutions.boolean),size(evolutions.pH,1),[]);
high_initial_pH_evolutions.alkalinity = reshape(evolutions.alkalinity(high_initial_pH_evolutions.boolean),size(evolutions.pH,1),[]);

for distribution_index = 1:size(high_initial_pH_evolutions.pH,1)
    high_initial_pH_evolutions.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],high_initial_pH_evolutions.pH(distribution_index,:));
    high_initial_pH_evolutions.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],high_initial_pH_evolutions.co2(distribution_index,:));
    high_initial_pH_evolutions.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,high_initial_pH_evolutions.saturation_state(distribution_index,:));
    high_initial_pH_evolutions.d11B_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:50,high_initial_pH_evolutions.d11B(distribution_index,:));
    high_initial_pH_evolutions.alkalinity_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:20000,high_initial_pH_evolutions.alkalinity(distribution_index,:)*1e6);
end
high_initial_pH_evolutions.pH_quantiles = high_initial_pH_evolutions.pH_distributions.quantile(quantiles);
high_initial_pH_evolutions.saturation_state_quantiles = high_initial_pH_evolutions.saturation_state_distributions.quantile(quantiles);
high_initial_pH_evolutions.d11B_quantiles = high_initial_pH_evolutions.d11B_distributions.quantile(quantiles);
high_initial_pH_evolutions.co2_quantiles = high_initial_pH_evolutions.co2_distributions.quantile(quantiles);

high_initial_pH_evolutions.pH_above_quantiles = high_initial_pH_evolutions.pH_distributions.quantile(above_quantiles);
high_initial_pH_evolutions.d11B_above_quantiles = high_initial_pH_evolutions.d11B_distributions.quantile(above_quantiles);
high_initial_pH_evolutions.co2_above_quantiles = high_initial_pH_evolutions.co2_distributions.quantile(above_quantiles);

%%
data_directory = "./../../Data/";
filename = "TJ_CO2_Evolutions_Metrics.csv";

writematrix([evolutions.pH_quantiles'],data_directory+filename,'WriteMode','overwrite');
writematrix([evolutions.saturation_state_quantiles'],data_directory+filename,'WriteMode','append');
writematrix([evolutions.co2_quantiles'],data_directory+filename,'WriteMode','append');

writematrix([high_initial_pH_evolutions.pH_quantiles'],data_directory+filename,'WriteMode','append');
writematrix([high_initial_pH_evolutions.saturation_state_quantiles'],data_directory+filename,'WriteMode','append');
writematrix([high_initial_pH_evolutions.co2_quantiles'],data_directory+filename,'WriteMode','append');
