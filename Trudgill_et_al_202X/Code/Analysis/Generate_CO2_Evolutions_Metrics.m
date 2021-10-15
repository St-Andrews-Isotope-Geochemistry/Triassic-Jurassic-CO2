clear all

%% Load data
boron_data = readtable("./../../Data/TJ_d11B.xlsx","Sheet","Delta_Temperature");
interpolation_ages = unique(sort([boron_data.age',linspace(min(boron_data.age),max(boron_data.age),80)]));

raw_evolutions = readmatrix("./../../Data/TJ_CO2_Evolutions.csv");
reshaped_evolutions = reshape(raw_evolutions,[100,13,numel(raw_evolutions)/(100*13)]);

evolutions.pH = squeeze(reshaped_evolutions(:,1,:));
evolutions.co2 = squeeze(reshaped_evolutions(:,2,:))*1e6;
evolutions.saturation_state = squeeze(reshaped_evolutions(:,3,:));
evolutions.dic = squeeze(reshaped_evolutions(:,4,:));
evolutions.alkalinity = squeeze(reshaped_evolutions(:,5,:));
evolutions.temperature = squeeze(reshaped_evolutions(:,6,:));
evolutions.d11B = squeeze(reshaped_evolutions(:,7,:));
evolutions.calcium = squeeze(reshaped_evolutions(:,8,:));
evolutions.magnesium = squeeze(reshaped_evolutions(:,9,:));
evolutions.epsilon = squeeze(reshaped_evolutions(:,10,:));
evolutions.d11B_sw = squeeze(reshaped_evolutions(:,11,:));

clear raw_evolutions reshaped_evolutions

is_nan = isnan(evolutions.pH);
evolutions.pH = reshape(evolutions.pH(~is_nan),numel(interpolation_ages),[]);
evolutions.co2 = reshape(evolutions.co2(~is_nan),numel(interpolation_ages),[]);
evolutions.saturation_state = reshape(evolutions.saturation_state(~is_nan),numel(interpolation_ages),[]);
evolutions.dic = reshape(evolutions.dic(~is_nan),numel(interpolation_ages),[]);
evolutions.alkalinity = reshape(evolutions.alkalinity(~is_nan),numel(interpolation_ages),[]);
evolutions.temperature = reshape(evolutions.temperature(~is_nan),numel(interpolation_ages),[]);
evolutions.d11B = reshape(evolutions.d11B(~is_nan),numel(interpolation_ages),[]);
evolutions.calcium = reshape(evolutions.calcium(~is_nan),numel(interpolation_ages),[]);
evolutions.magnesium = reshape(evolutions.magnesium(~is_nan),numel(interpolation_ages),[]);
evolutions.epsilon = reshape(evolutions.epsilon(~is_nan),numel(interpolation_ages),[]);
evolutions.d11B_sw = reshape(evolutions.d11B_sw(~is_nan),numel(interpolation_ages),[]);

quantile_level = 5; % As percent
quantiles = [(quantile_level/100)/2,0.5,1-((quantile_level/100)/2)];
above_quantiles = [0.05,0.5,0.95];

%% Get the initial subsample
for distribution_index = 1:size(evolutions.pH,1)
    evolutions.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],evolutions.pH(distribution_index,:));
    evolutions.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],evolutions.co2(distribution_index,:));
    evolutions.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,evolutions.saturation_state(distribution_index,:));
    evolutions.d11B_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:50,evolutions.d11B(distribution_index,:));
    evolutions.alkalinity_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:20000,evolutions.alkalinity(distribution_index,:)*1e6);
end
evolutions.pH_quantiles = evolutions.pH_distributions.quantile(quantiles);
evolutions.d11B_quantiles = evolutions.d11B_distributions.quantile(quantiles);
evolutions.co2_quantiles = evolutions.co2_distributions.quantile(quantiles);

evolutions.pH_above_quantiles = evolutions.pH_distributions.quantile(above_quantiles);
evolutions.d11B_above_quantiles = evolutions.d11B_distributions.quantile(above_quantiles);
evolutions.co2_above_quantiles = evolutions.co2_distributions.quantile(above_quantiles);

%% High pH subsample
high_initial_co2.boolean = repmat(evolutions.co2(end,:)>=800,size(evolutions.pH,1),1);
high_initial_co2.pH = reshape(evolutions.pH(high_initial_co2.boolean),100,[]);
high_initial_co2.co2 =  reshape(evolutions.co2(high_initial_co2.boolean),100,[]);
high_initial_co2.saturation_state = reshape(evolutions.saturation_state(high_initial_co2.boolean),100,[]);
high_initial_co2.d11B = reshape(evolutions.d11B(high_initial_co2.boolean),100,[]);

for distribution_index = 1:size(high_initial_co2.pH,1)
    high_initial_co2.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],high_initial_co2.pH(distribution_index,:));
    high_initial_co2.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],high_initial_co2.co2(distribution_index,:));
    high_initial_co2.d11B_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],high_initial_co2.d11B(distribution_index,:));
    high_initial_co2.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],high_initial_co2.saturation_state(distribution_index,:));
end

high_initial_co2.pH_quantiles = high_initial_co2.pH_distributions.quantile(quantiles);
high_initial_co2.d11B_quantiles = high_initial_co2.d11B_distributions.quantile(quantiles);
high_initial_co2.co2_quantiles = high_initial_co2.co2_distributions.quantile(quantiles);
high_initial_co2.saturation_state_quantiles = high_initial_co2.saturation_state_distributions.quantile(quantiles);

high_initial_co2.pH_above_quantiles = high_initial_co2.pH_distributions.quantile(above_quantiles);
high_initial_co2.d11B_above_quantiles = high_initial_co2.d11B_distributions.quantile(above_quantiles);
high_initial_co2.co2_above_quantiles = high_initial_co2.co2_distributions.quantile(above_quantiles);

%% Medium CO2 subsample
medium_initial_co2.boolean = repmat(evolutions.co2(end,:)>=600 & evolutions.co2(end,:)<=2000,size(evolutions.pH,1),1);
medium_initial_co2.pH = reshape(evolutions.pH(medium_initial_co2.boolean),100,[]);
medium_initial_co2.co2 =  reshape(evolutions.co2(medium_initial_co2.boolean),100,[]);
medium_initial_co2.saturation_state = reshape(evolutions.saturation_state(medium_initial_co2.boolean),100,[]);
medium_initial_co2.d11B = reshape(evolutions.d11B(medium_initial_co2.boolean),100,[]);

for distribution_index = 1:size(medium_initial_co2.pH,1)
    medium_initial_co2.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],medium_initial_co2.pH(distribution_index,:));
    medium_initial_co2.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],medium_initial_co2.co2(distribution_index,:));
    medium_initial_co2.d11B_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],medium_initial_co2.d11B(distribution_index,:));
    medium_initial_co2.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],medium_initial_co2.saturation_state(distribution_index,:));
end
medium_initial_co2.pH_quantiles = medium_initial_co2.pH_distributions.quantile(quantiles);
medium_initial_co2.d11B_quantiles = medium_initial_co2.d11B_distributions.quantile(quantiles);
medium_initial_co2.co2_quantiles = medium_initial_co2.co2_distributions.quantile(quantiles);
medium_initial_co2.saturation_state_quantiles = medium_initial_co2.saturation_state_distributions.quantile(quantiles);

medium_initial_co2.pH_above_quantiles = medium_initial_co2.pH_distributions.quantile(above_quantiles);
medium_initial_co2.d11B_above_quantiles = medium_initial_co2.d11B_distributions.quantile(above_quantiles);
medium_initial_co2.co2_above_quantiles = medium_initial_co2.co2_distributions.quantile(above_quantiles);

%% Low pH subsample
low_initial_co2.boolean = repmat(evolutions.co2(end,:)<800,size(evolutions.pH,1),1);
low_initial_co2.pH = reshape(evolutions.pH(low_initial_co2.boolean),100,[]);
low_initial_co2.co2 =  reshape(evolutions.co2(low_initial_co2.boolean),100,[]);
low_initial_co2.saturation_state = reshape(evolutions.saturation_state(low_initial_co2.boolean),100,[]);
low_initial_co2.d11B = reshape(evolutions.d11B(low_initial_co2.boolean),100,[]);

for distribution_index = 1:size(low_initial_co2.pH,1)
    low_initial_co2.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],low_initial_co2.pH(distribution_index,:));
    low_initial_co2.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],low_initial_co2.co2(distribution_index,:));
    low_initial_co2.d11B_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],low_initial_co2.d11B(distribution_index,:));
    low_initial_co2.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],low_initial_co2.saturation_state(distribution_index,:));
    low_initial_co2.d11B_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],low_initial_co2.d11B(distribution_index,:));
end

low_initial_co2.pH_quantiles = low_initial_co2.pH_distributions.quantile(quantiles);
low_initial_co2.d11B_quantiles = low_initial_co2.d11B_distributions.quantile(quantiles);
low_initial_co2.co2_quantiles = low_initial_co2.co2_distributions.quantile(quantiles);
low_initial_co2.saturation_state_quantiles = low_initial_co2.saturation_state_distributions.quantile(quantiles);

low_initial_co2.pH_above_quantiles = low_initial_co2.pH_distributions.quantile(above_quantiles);
low_initial_co2.d11B_above_quantiles = low_initial_co2.d11B_distributions.quantile(above_quantiles);
low_initial_co2.co2_above_quantiles = low_initial_co2.co2_distributions.quantile(above_quantiles);

%% Temperature
preperturbation.boolean = (interpolation_ages>=boron_data.age(9))';
preperturbation.temperature_evolutions = evolutions.temperature(preperturbation.boolean,:);
preperturbation.temperature_mean = mean(preperturbation.temperature_evolutions);
preperturbation.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(20,preperturbation.temperature_mean);

perturbation.boolean = (interpolation_ages==boron_data.age(19) | interpolation_ages==boron_data.age(20))';
perturbation.temperature_evolutions = evolutions.temperature(perturbation.boolean,:);
perturbation.temperature_mean = mean(perturbation.temperature_evolutions);
perturbation.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(20,perturbation.temperature_mean);

temperature_change = perturbation.temperature_mean-preperturbation.temperature_mean;
temperature_change_distribution = Geochemistry_Helpers.Distribution.fromSamples([],temperature_change);

%%
data_directory = "./../../Data/";
filename = "TJ_CO2_Evolutions_Metrics.csv";

writematrix([low_initial_co2.pH_quantiles';medium_initial_co2.pH_quantiles';high_initial_co2.pH_quantiles'],data_directory+filename,'WriteMode','overwrite');
writematrix([low_initial_co2.saturation_state_quantiles';medium_initial_co2.saturation_state_quantiles';high_initial_co2.saturation_state_quantiles'],data_directory+filename,'WriteMode','append');
writematrix([low_initial_co2.co2_quantiles';medium_initial_co2.co2_quantiles';high_initial_co2.co2_quantiles'],data_directory+filename,'WriteMode','append');


