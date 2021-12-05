clear

%% Load data
d18O_d13C = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Delta_Temperature");
d18O_d13C_averaged = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Delta_Temperature");

boron_data = readtable("./../../Data/TJ_d11B.xlsx","Sheet","Delta_Temperature");
evolutions = getValidSamples(getShapedEvolutions("./../../Data/TJ_CO2_Evolutions.csv"));

interpolation_ages = unique(sort([boron_data.age',linspace(min(boron_data.age),max(boron_data.age),80)]));
evolutions.age = repmat(interpolation_ages',1,size(evolutions.pH,2));

[preperturbation_evolutions,perturbation_evolutions] = splitPerturbation(evolutions,boron_data.age(9));

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