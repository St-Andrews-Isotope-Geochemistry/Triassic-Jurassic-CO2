clear

%% Load data
d18O_d13C = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab");
d18O_d13C_averaged = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Averaged");

boron_data = readtable("./../../Data/TJ_d11B_pH.xlsx");
boron_data.age = boron_data.absolute_age;

raw_evolutions = readmatrix("./../../Data/TJ_CO2_Evolutions.csv");
reshaped_evolutions = reshape(raw_evolutions,[22,11,100000]);

evolutions.pH = squeeze(reshaped_evolutions(:,1,:));
evolutions.co2 = squeeze(reshaped_evolutions(:,2,:));
evolutions.saturation_state = squeeze(reshaped_evolutions(:,3,:));
evolutions.dic = squeeze(reshaped_evolutions(:,4,:));
evolutions.alkalinity = squeeze(reshaped_evolutions(:,5,:));
evolutions.temperature = squeeze(reshaped_evolutions(:,6,:));
evolutions.d11B = squeeze(reshaped_evolutions(:,7,:));
evolutions.calcium = squeeze(reshaped_evolutions(:,8,:));
evolutions.magnesium = squeeze(reshaped_evolutions(:,9,:));
evolutions.epsilon = squeeze(reshaped_evolutions(:,10,:));
evolutions.d11B_sw = squeeze(reshaped_evolutions(:,1,:));

clear raw_evolutions reshaped_evolutions

%% Get the initial subsample
evolutions.subsample_boolean = repmat(evolutions.saturation_state(1,:)>=5 & evolutions.saturation_state(1,:)<=10.7 & evolutions.co2(1,:)>=400 & evolutions.co2(1,:)<=5000 & all(evolutions.co2<10000) & all(evolutions.co2>0) & all(evolutions.saturation_state<12),size(evolutions.pH,1),1);
evolutions.pH_subsample = reshape(evolutions.pH(evolutions.subsample_boolean),22,[]);
evolutions.co2_subsample =  reshape(evolutions.co2(evolutions.subsample_boolean),22,[]);
evolutions.saturation_state_subsample = reshape(evolutions.saturation_state(evolutions.subsample_boolean),22,[]);

%% High pH subsample
high_initial_co2.boolean = repmat(evolutions.co2_subsample(1,:)>=2000,size(evolutions.pH_subsample,1),1);
high_initial_co2.pH = reshape(evolutions.pH_subsample(high_initial_co2.boolean),22,[]);
high_initial_co2.co2 =  reshape(evolutions.co2_subsample(high_initial_co2.boolean),22,[]);
high_initial_co2.saturation_state = reshape(evolutions.saturation_state_subsample(high_initial_co2.boolean),22,[]);

for distribution_index = 1:size(high_initial_co2.pH,1)
    high_initial_co2.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,high_initial_co2.pH(distribution_index,:));
    high_initial_co2.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,high_initial_co2.co2(distribution_index,:));
    high_initial_co2.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,high_initial_co2.saturation_state(distribution_index,:));
end

high_initial_co2.pH_quantiles = [high_initial_co2.pH_distributions.quantile(0.025),high_initial_co2.pH_distributions.quantile(0.5),high_initial_co2.pH_distributions.quantile(0.975)];
high_initial_co2.saturation_state_quantiles = [high_initial_co2.saturation_state_distributions.quantile(0.025),high_initial_co2.saturation_state_distributions.quantile(0.5),high_initial_co2.saturation_state_distributions.quantile(0.975)];
high_initial_co2.co2_quantiles = [high_initial_co2.co2_distributions.quantile(0.025),high_initial_co2.co2_distributions.quantile(0.5),high_initial_co2.co2_distributions.quantile(0.975)];

%% Medium CO2 subsample
medium_initial_co2.boolean = repmat(evolutions.co2_subsample(1,:)>=600 & evolutions.co2_subsample(1,:)<=2000,size(evolutions.pH_subsample,1),1);
medium_initial_co2.pH = reshape(evolutions.pH_subsample(medium_initial_co2.boolean),22,[]);
medium_initial_co2.co2 =  reshape(evolutions.co2_subsample(medium_initial_co2.boolean),22,[]);
medium_initial_co2.saturation_state = reshape(evolutions.saturation_state_subsample(medium_initial_co2.boolean),22,[]);

for distribution_index = 1:size(medium_initial_co2.pH,1)
    medium_initial_co2.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,medium_initial_co2.pH(distribution_index,:));
    medium_initial_co2.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,medium_initial_co2.co2(distribution_index,:));
    medium_initial_co2.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,medium_initial_co2.saturation_state(distribution_index,:));
end

medium_initial_co2.pH_quantiles = [medium_initial_co2.pH_distributions.quantile(0.025),medium_initial_co2.pH_distributions.quantile(0.5),medium_initial_co2.pH_distributions.quantile(0.975)];
medium_initial_co2.saturation_state_quantiles = [medium_initial_co2.saturation_state_distributions.quantile(0.025),medium_initial_co2.saturation_state_distributions.quantile(0.5),medium_initial_co2.saturation_state_distributions.quantile(0.975)];
medium_initial_co2.co2_quantiles = [medium_initial_co2.co2_distributions.quantile(0.025),medium_initial_co2.co2_distributions.quantile(0.5),medium_initial_co2.co2_distributions.quantile(0.975)];

%% Low pH subsample
low_initial_co2.boolean = repmat(evolutions.co2_subsample(1,:)>400 & evolutions.co2_subsample(1,:)<600,size(evolutions.pH_subsample,1),1);
low_initial_co2.pH = reshape(evolutions.pH_subsample(low_initial_co2.boolean),22,[]);
low_initial_co2.co2 =  reshape(evolutions.co2_subsample(low_initial_co2.boolean),22,[]);
low_initial_co2.saturation_state = reshape(evolutions.saturation_state_subsample(low_initial_co2.boolean),22,[]);

for distribution_index = 1:size(low_initial_co2.pH,1)
    low_initial_co2.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,low_initial_co2.pH(distribution_index,:));
    low_initial_co2.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,low_initial_co2.co2(distribution_index,:));
    low_initial_co2.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,low_initial_co2.saturation_state(distribution_index,:));
end

low_initial_co2.pH_quantiles = [low_initial_co2.pH_distributions.quantile(0.025),low_initial_co2.pH_distributions.quantile(0.5),low_initial_co2.pH_distributions.quantile(0.975)];
low_initial_co2.saturation_state_quantiles = [low_initial_co2.saturation_state_distributions.quantile(0.025),low_initial_co2.saturation_state_distributions.quantile(0.5),low_initial_co2.saturation_state_distributions.quantile(0.975)];
low_initial_co2.co2_quantiles = [low_initial_co2.co2_distributions.quantile(0.025),low_initial_co2.co2_distributions.quantile(0.5),low_initial_co2.co2_distributions.quantile(0.975)];


%%
data_directory = "./../../Data/";
filename = "TJ_CO2_Evolutions_Metrics.csv";


writematrix([low_initial_co2.pH_quantiles';medium_initial_co2.pH_quantiles';high_initial_co2.pH_quantiles'],data_directory+filename,'WriteMode','overwrite');
writematrix([low_initial_co2.saturation_state_quantiles';medium_initial_co2.saturation_state_quantiles';high_initial_co2.saturation_state_quantiles'],data_directory+filename,'WriteMode','append');
writematrix([low_initial_co2.co2_quantiles';medium_initial_co2.co2_quantiles';high_initial_co2.co2_quantiles'],data_directory+filename,'WriteMode','append');


