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
evolutions.subsample_boolean = repmat(evolutions.saturation_state(1,:)>=5 & evolutions.saturation_state(1,:)<=10.7 & evolutions.co2(1,:)>=400 & evolutions.co2(1,:)<=5000 & all(evolutions.co2>0) & all(evolutions.saturation_state<12),size(evolutions.pH,1),1);
evolutions.pH_subsample = reshape(evolutions.pH(evolutions.subsample_boolean),22,[]);
evolutions.co2_subsample =  reshape(evolutions.co2(evolutions.subsample_boolean),22,[]);
evolutions.saturation_state_subsample = reshape(evolutions.saturation_state(evolutions.subsample_boolean),22,[]);

%% Plot for the distribution of CO2 change
co2_change = max(evolutions.co2_subsample(10:end,:))-mean(evolutions.co2_subsample (1:9,:));
co2_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(-1000:100:100000,co2_change).normalise();

co2_change_quantiles = [co2_change_distribution.quantile(0.05),co2_change_distribution.quantile(0.5),co2_change_distribution.quantile(0.95)];

figure(1);
clf
hold on
co2_change_distribution.plot();
plot([co2_change_quantiles;co2_change_quantiles],[0;1],'k');

ylim([0,max(co2_change_distribution.probabilities)]);
xlabel("CO_2 (ppm)");
ylabel("Probability");

current_axis = gca;
current_axis.XAxis.Exponent = 0;