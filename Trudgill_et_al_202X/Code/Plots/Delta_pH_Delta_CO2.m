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

%%
evolutions.delta_pH = diff(evolutions.pH);
evolutions.delta_co2 = diff(evolutions.co2);

%%
delta_pH_distributions = Geochemistry_Helpers.Distribution().create(numel(boron_data.age)-1);
delta_co2_distributions = Geochemistry_Helpers.Distribution().create(numel(boron_data.age)-1);
for age_index = 1:numel(boron_data.age)-1
    delta_pH_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(-1:0.01:1,evolutions.delta_pH(age_index,:)).normalise();
    delta_co2_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(-10000:100:10000,evolutions.delta_co2(age_index,:)).normalise();
end

%
delta_pH_median = delta_pH_distributions.median();
delta_pH_uncertainties = [delta_pH_distributions.quantile(0.025),delta_pH_distributions.quantile(0.975)];

delta_co2_median = delta_co2_distributions.median();
delta_co2_uncertainties = [delta_co2_distributions.quantile(0.025),delta_co2_distributions.quantile(0.975)];

%% Get the initial subsample
evolutions_low_co2.subsample_boolean = repmat(evolutions.co2(1,:)<=2000 & all(evolutions.co2>0) & all(evolutions.saturation_state<12),size(evolutions.pH,1),1);
evolutions_low_co2.pH = reshape(evolutions.pH(evolutions_low_co2.subsample_boolean),22,[]);
evolutions_low_co2.co2 =  reshape(evolutions.co2(evolutions_low_co2.subsample_boolean),22,[]);
% evolutions.saturation_state_subsample = reshape(evolutions.saturation_state(evolutions.subsample_boolean),22,[]);

evolutions_low_co2.delta_pH = diff(evolutions_low_co2.pH);
evolutions_low_co2.delta_co2 = diff(evolutions_low_co2.co2);

evolutions_low_co2.delta_pH_distributions = Geochemistry_Helpers.Distribution().create(numel(boron_data.age)-1);
evolutions_low_co2.delta_co2_distributions = Geochemistry_Helpers.Distribution().create(numel(boron_data.age)-1);
for age_index = 1:numel(boron_data.age)-1
    evolutions_low_co2.delta_pH_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(-1:0.01:1,evolutions_low_co2.delta_pH(age_index,:)).normalise();
    evolutions_low_co2.delta_co2_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(-10000:100:10000,evolutions_low_co2.delta_co2(age_index,:)).normalise();
end

%
evolutions_low_co2.delta_pH_median = evolutions_low_co2.delta_pH_distributions.median();
evolutions_low_co2.delta_pH_uncertainties = [evolutions_low_co2.delta_pH_distributions.quantile(0.025),evolutions_low_co2.delta_pH_distributions.quantile(0.975)];

evolutions_low_co2.delta_co2_median = evolutions_low_co2.delta_co2_distributions.median();
evolutions_low_co2.delta_co2_uncertainties = [evolutions_low_co2.delta_co2_distributions.quantile(0.025),evolutions_low_co2.delta_co2_distributions.quantile(0.975)];

%% 
figure(1);
clf
subplot(2,1,1);
hold on
patch([boron_data.age(1:end-1);flipud(boron_data.age(1:end-1))],[delta_pH_uncertainties(:,1);flipud(delta_pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'EdgeColor','None','FaceAlpha',0.5);
% patch([boron_data.age(1:end-1);flipud(boron_data.age(1:end-1))],[evolutions_low_co2.delta_pH_uncertainties(:,1);flipud(evolutions_low_co2.delta_pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("ForestGreen").rgb,'EdgeColor','None','FaceAlpha',0.5);

plot(boron_data.age(1:end-1),delta_pH_median);

set(gca,'XDir','Reverse','XTick',[]);
ylabel("\DeltapH");

subplot(2,1,2);
hold on
patch([boron_data.age(1:end-1);flipud(boron_data.age(1:end-1))],[delta_co2_uncertainties(:,1);flipud(delta_co2_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("FireBrick").rgb,'EdgeColor','None','FaceAlpha',0.5);
plot(boron_data.age(1:end-1),delta_co2_median,'Color',Geochemistry_Helpers.Colour.Colour("Crimson").rgb);

current_axis = gca;
current_axis.YAxis.Exponent = 0;
set(gca,'XDir','Reverse');
xlabel("Age (Ma)");
ylabel("\DeltaCO_2");

exportgraphics(gcf,"./../../Figures/Delta_pH_Delta_CO2.png");
