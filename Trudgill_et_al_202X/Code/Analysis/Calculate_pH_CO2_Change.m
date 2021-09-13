clear

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

preperturbation_boundary = boron_data.age(9);

%% Get the initial subsample
for distribution_index = 1:size(evolutions.pH,1)
    evolutions.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(4:0.01:11,evolutions.pH(distribution_index,:));
    evolutions.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,evolutions.co2(distribution_index,:));
    evolutions.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,evolutions.saturation_state(distribution_index,:));
    evolutions.d11B_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:50,evolutions.d11B(distribution_index,:));
    evolutions.alkalinity_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:20000,evolutions.alkalinity(distribution_index,:)*1e6);
end
evolutions.pH_quantiles = evolutions.pH_distributions.quantile(quantiles);
evolutions.d11B_quantiles = evolutions.d11B_distributions.quantile(quantiles);

%%
preperturbation.pH = evolutions.pH(interpolation_ages>=preperturbation_boundary,:);
preperturbation.co2 = evolutions.co2(interpolation_ages>=preperturbation_boundary,:);
preperturbation.pH_mean = mean(preperturbation.pH,1,'omitnan');
perperturbation.pH_distribution = Geochemistry_Helpers.Distribution.fromSamples(4:0.05:11,preperturbation.pH_mean);
preperturbation.co2_mean = mean(preperturbation.co2,1,'omitnan');
perperturbation.co2_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:100:50000,preperturbation.co2_mean);

perturbation.pH = evolutions.pH(interpolation_ages<preperturbation_boundary,:);
perturbation.co2 = evolutions.co2(interpolation_ages<preperturbation_boundary,:);
perturbation.pH_minimum = min(perturbation.pH,[],1,'omitnan');
perturbation.pH_distribution = Geochemistry_Helpers.Distribution.fromSamples(4:0.05:11,perturbation.pH_minimum);
perturbation.co2_maximum = max(perturbation.co2,[],1,'omitnan');
perturbation.co2_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:100:50000,perturbation.co2_maximum);

%%
pH_change = perturbation.pH_minimum-preperturbation.pH_mean;
co2_change = perturbation.co2_maximum-preperturbation.co2_mean;

pH_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(-3:0.05:3,pH_change);
co2_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(100,co2_change);

%%
figure(1);
clf
hold on
perperturbation.pH_distribution.plot();
perturbation.pH_distribution.plot();

figure(2);
clf
hold on
perperturbation.co2_distribution.plot();
perturbation.co2_distribution.plot();

figure(3);
clf
hold on
pH_change_distribution.plot();

figure(4);
clf
hold on
co2_change_distribution.plot();
