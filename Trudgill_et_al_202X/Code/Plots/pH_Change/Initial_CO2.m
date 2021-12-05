clear

%% Load data
d18O_d13C = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Delta_Temperature");
d18O_d13C_averaged = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Averaged");

boron_data = readtable("./../../Data/TJ_d11B.xlsx","Sheet","Delta_Temperature");
evolutions = getShapedEvolutions("./../../Data/TJ_CO2_Evolutions.csv");

%% Get the initial subsample
evolutions_low_co2.subsample_boolean = repmat(evolutions.co2(1,:)<=2000 & all(evolutions.co2>0) & all(evolutions.saturation_state<12),size(evolutions.pH,1),1);
evolutions_low_co2.pH = reshape(evolutions.pH(evolutions_low_co2.subsample_boolean),100,[]);
evolutions_low_co2.co2 =  reshape(evolutions.co2(evolutions_low_co2.subsample_boolean),100,[]);

evolutions_low_co2.pH_distributions = Geochemistry_Helpers.Distribution().create(numel(boron_data.age)-1);
evolutions_low_co2.co2_distributions = Geochemistry_Helpers.Distribution().create(numel(boron_data.age)-1);
for age_index = 1:numel(boron_data.age)
    evolutions_low_co2.pH_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(5:0.01:10,evolutions_low_co2.pH(age_index,:)).normalise();
    evolutions_low_co2.co2_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,evolutions_low_co2.co2(age_index,:)).normalise();
end

evolutions_low_co2.pH_median = evolutions_low_co2.pH_distributions.median();
evolutions_low_co2.pH_uncertainties = [evolutions_low_co2.pH_distributions.quantile(0.025),evolutions_low_co2.pH_distributions.quantile(0.975)];

evolutions_low_co2.co2_median = evolutions_low_co2.co2_distributions.median();
evolutions_low_co2.co2_uncertainties = [evolutions_low_co2.co2_distributions.quantile(0.025),evolutions_low_co2.co2_distributions.quantile(0.975)];

%
evolutions_high_co2.subsample_boolean = repmat(evolutions.co2(1,:)>2000 & all(evolutions.co2>0) & all(evolutions.saturation_state<12),size(evolutions.pH,1),1);
evolutions_high_co2.pH = reshape(evolutions.pH(evolutions_high_co2.subsample_boolean),100,[]);
evolutions_high_co2.co2 =  reshape(evolutions.co2(evolutions_high_co2.subsample_boolean),100,[]);

evolutions_high_co2.pH_distributions = Geochemistry_Helpers.Distribution().create(numel(boron_data.age)-1);
evolutions_high_co2.co2_distributions = Geochemistry_Helpers.Distribution().create(numel(boron_data.age)-1);
for age_index = 1:numel(boron_data.age)
    evolutions_high_co2.pH_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(5:0.01:10,evolutions_high_co2.pH(age_index,:)).normalise();
    evolutions_high_co2.co2_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,evolutions_high_co2.co2(age_index,:)).normalise();
end

evolutions_high_co2.pH_median = evolutions_high_co2.pH_distributions.median();
evolutions_high_co2.pH_uncertainties = [evolutions_high_co2.pH_distributions.quantile(0.025),evolutions_high_co2.pH_distributions.quantile(0.975)];

evolutions_high_co2.co2_median = evolutions_high_co2.co2_distributions.median();
evolutions_high_co2.co2_uncertainties = [evolutions_high_co2.co2_distributions.quantile(0.025),evolutions_high_co2.co2_distributions.quantile(0.975)];

%%
figure(1);
clf

subplot_handles{1} = subplot(2,1,1);
hold on
plot_handles{1} = patch([boron_data.age;flipud(boron_data.age)],[evolutions_low_co2.pH_uncertainties(:,1);flipud(evolutions_low_co2.pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot_handles{2} = patch([boron_data.age;flipud(boron_data.age)],[evolutions_high_co2.pH_uncertainties(:,1);flipud(evolutions_high_co2.pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("Crimson").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("Crimson").darken(1.2).rgb);

set(gca,'XDir','Reverse','XTick',[]);
ylabel("pH");

subplot_handles{2} = subplot(2,1,2);
plot_handles{1} = patch([boron_data.age;flipud(boron_data.age)],[evolutions_low_co2.co2_uncertainties(:,1);flipud(evolutions_low_co2.co2_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot_handles{2} = patch([boron_data.age;flipud(boron_data.age)],[evolutions_high_co2.co2_uncertainties(:,1);flipud(evolutions_high_co2.co2_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("Crimson").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("Crimson").darken(1.2).rgb);

legend([plot_handles{:}],["Initial CO_2 < 2000ppm","Initial CO_2 > 2000ppm"],'Location','NorthWest','Box','Off');

set(gca,'XDir','Reverse');
xlabel("Age (Ma)");
ylabel("CO_2 (ppm)");
ylim([0,11000]);

exportgraphics(gcf,"./../../Figures/Initial_CO2.png");

