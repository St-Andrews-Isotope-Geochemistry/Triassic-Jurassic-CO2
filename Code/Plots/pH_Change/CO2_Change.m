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
preperturbation_evolutions.co2_mean = mean(preperturbation_evolutions.co2);

perturbation_evolutions.d11B_minimum_boolean = perturbation_evolutions.d11B == min(perturbation_evolutions.d11B);
% evolutions_d11B_minimum.subsample_boolean = d11B_minimum;
perturbation_evolutions.d11B_minimum_pH = perturbation_evolutions.pH(perturbation_evolutions.d11B_minimum_boolean)';
perturbation_evolutions.d11B_minimum_co2 =  perturbation_evolutions.co2(perturbation_evolutions.d11B_minimum_boolean)';

co2_change = (perturbation_evolutions.d11B_minimum_co2-preperturbation_evolutions.co2_mean);
co2_change_distribution = Geochemistry_Helpers.Distribution().fromSamples(-10000:100:100000,co2_change);

co2_fractional_change = co2_change./preperturbation_evolutions.co2_mean;
co2_doublings = log2(co2_fractional_change);
% co2_doublings = co2_doublings(imag(co2_doublings)==0);

co2_doublings_distribution = Geochemistry_Helpers.Distribution().fromSamples(0:0.1:6,co2_doublings);

%%
figure(1);
clf
co2_change_distribution.plot();

set(gca,'XTick',[-2:2:10]*10^4,'XTickLabel',num2str([-2:2:10]'*10^4));
xlabel("\DeltaCO_2 (ppm)");
ylabel("Probability");

figure(2);
clf
co2_doublings_distribution.plot();

xlabel("\DeltaCO_2 (doublings)");
ylabel("Probability");