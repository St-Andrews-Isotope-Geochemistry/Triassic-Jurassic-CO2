% TJ_CO2
% Using: CO2 between 500ppm and 5000ppm
%        Omega between 5 and 10.7
% Gives: Viable initial alkalinity distribution
%        Viable initial pH distribution
%
% Using: Initial pH distribution
%        Measured d11B
%        Epsilon + Kb
% Gives: d11B_sw distribution
%
% Using: Measured d11B range
% Gives: d11B_sw distribution
%
% Combine both d11B_sw distributions
%
% Using: Assumptions on rate of change of alkalinity
% Gives: Possible alkalinity evolutions
%        
% Using: d11B_sw distribution
%        Measured d11B
%        Possible alkalinity evolutions
% Gives: pH distribution evolution
%        CO2 distribution evolution

clear
%%
tic
% Choose the number of statistical samples
number_of_samples = 10000;

% Load in the data
temperature_data = readtable("./../Data/TJ_Temperature.xlsx");
boron_data = readtable("./../Data/TJ_d11B_pH.xlsx");

boron_data.age = boron_data.absolute_age;

% Segregate into before + during/after perturbation
background_data = boron_data(1:9,:);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4 = min(perturbation_data.d11B,[],'omitnan');

temperature_change = [0;diff(boron_data.hansen_temperature)];

%% Distributions and Samplers
% Generate distributions for each variable to represent uncertainties
initial_co2_distribution = Geochemistry_Helpers.Distribution(0:100e-6:10000e-6,"Flat",[400e-6,5000e-6]).normalise();
initial_omega_distribution = Geochemistry_Helpers.Distribution(0:0.1:12,"Flat",[5,10.7]).normalise();
initial_temperature_distribution = Geochemistry_Helpers.Distribution(-10:5:60,"Gaussian",[14.8,3.3]).normalise();
temperature_change_distribution = Geochemistry_Helpers.Distribution(0:0.1:10,"Gaussian",[4.2,1.1]).normalise();
ca_distribution = Geochemistry_Helpers.Distribution(0:0.1:20,"Flat",[8,17]).normalise();
mg_distribution = Geochemistry_Helpers.Distribution(20:0.1:61,"Flat",[28,52]).normalise();
epsilon_distribution = Geochemistry_Helpers.Distribution(27:0.01:31,"Gaussian",[27.2,0.6]).normalise();
d11B_distributions = Geochemistry_Helpers.Distribution().create([height(boron_data),1]);

for d11B_distribution_index = 1:numel(d11B_distributions)
    d11B_distributions(d11B_distribution_index) = Geochemistry_Helpers.Distribution(0:0.1:40,"Gaussian",[boron_data.d11B(d11B_distribution_index),boron_data.uncertainty(d11B_distribution_index)/2]).normalise();
end

% Use a sampling technique for those distributions
initial_co2_sampler = Geochemistry_Helpers.Sampler(initial_co2_distribution,"latin_hypercube");
initial_omega_sampler = Geochemistry_Helpers.Sampler(initial_omega_distribution,"latin_hypercube");
initial_temperature_sampler = Geochemistry_Helpers.Sampler(initial_temperature_distribution,"latin_hypercube");
temperature_change_sampler = Geochemistry_Helpers.Sampler(temperature_change_distribution,"latin_hypercube");
ca_sampler = Geochemistry_Helpers.Sampler(ca_distribution,"latin_hypercube");
mg_sampler = Geochemistry_Helpers.Sampler(mg_distribution,"latin_hypercube");
epsilon_sampler = Geochemistry_Helpers.Sampler(epsilon_distribution,"latin_hypercube");

d11B_samplers = Geochemistry_Helpers.Sampler().create([height(boron_data),1]);
for d11B_sampler_index = 1:numel(d11B_distributions)
    d11B_samplers(d11B_sampler_index) = Geochemistry_Helpers.Sampler(d11B_distributions(d11B_sampler_index),"latin_hypercube");
end

% Get the samples
initial_co2_sampler.getSamples(number_of_samples).shuffle();
initial_omega_sampler.getSamples(number_of_samples).shuffle();
initial_temperature_sampler.getSamples(number_of_samples).shuffle();
temperature_change_sampler.getSamples(number_of_samples).shuffle();
ca_sampler.getSamples(number_of_samples).shuffle();
mg_sampler.getSamples(number_of_samples).shuffle();
epsilon_sampler.getSamples(number_of_samples).shuffle();

for d11B_sampler_index = 1:numel(d11B_samplers)
    d11B_samplers(d11B_sampler_index).getSamples(number_of_samples).shuffle();
end
d11B_evolutions = d11B_samplers.collate("samples");

initial_d11B_CO2 = BuCC.d11BCO2().create(number_of_samples);

% Species calibration
initial_d11B_CO2.collate("species_calibration").collate("d11B_measured").assignToEach("value",d11B_samplers(1).samples);

% Boron
initial_d11B_CO2.collate("boron").collate("pH").assignToAll("pValue",NaN);
initial_d11B_CO2.collate("boron").collate("d11B_sw").assignToAll("value",NaN);
initial_d11B_CO2.collate("boron").assignToEach("epsilon",epsilon_sampler.samples);

% Carbonate chemistry
% Main parameters
initial_d11B_CO2.collate("carbonate_chemistry").collate("atmospheric_co2").assignToEach("partial_pressure",initial_co2_sampler.samples);
initial_d11B_CO2.collate("carbonate_chemistry").assignToEach("saturation_state",initial_omega_sampler.samples);
initial_d11B_CO2.collate("carbonate_chemistry").assignToAll("units"," mol/kg");

% Ancillary
% Variable
initial_d11B_CO2.collate("carbonate_chemistry").assignToEach("temperature",initial_temperature_sampler.samples);
initial_d11B_CO2.collate("carbonate_chemistry").assignToEach("calcium",ca_sampler.samples);
initial_d11B_CO2.collate("carbonate_chemistry").assignToEach("magnesium",mg_sampler.samples);

% Constant
initial_d11B_CO2.collate("carbonate_chemistry").assignToAll("salinity",35);
initial_d11B_CO2.collate("carbonate_chemistry").assignToAll("oceanic_pressure",0);
initial_d11B_CO2.collate("carbonate_chemistry").assignToAll("atmospheric_pressure",1);

% Create a MyAMI object
myami = MyAMI.MyAMI("Precalculated",true);
initial_d11B_CO2.collate("carbonate_chemistry").collate("equilibrium_coefficients").assignToAll("MyAMI",myami);

% Do carbonate chemistry calculations
initial_d11B_CO2.calculate();

% Collate the results of those calculations back to normal arrays
initial_pH = initial_d11B_CO2.collate("carbonate_chemistry").collate("pH").collate("pValue");
initial_alkalinity = initial_d11B_CO2.collate("carbonate_chemistry").collate("alkalinity")*1e6; %(converted units from mol/kg to umol/kg)

% Create a distribution from those sample arrays
initial_pH_distribution = Geochemistry_Helpers.Distribution.fromSamples((6:0.01:9)',initial_pH).normalise();
initial_alkalinity_distribution = Geochemistry_Helpers.Distribution.fromSamples((0:10:16000)',initial_alkalinity).normalise();

% Smooth the distribution of alkalinity to minimise noise (from statistical
% randomness)
initial_alkalinity_smoothed = smooth(initial_alkalinity_distribution.probabilities,sqrt(number_of_samples));
% Create a new distribution from the smoothed data and a sampler to get
% plausible alkalinities
initial_alkalinity_smoothed_distribution = Geochemistry_Helpers.Distribution(initial_alkalinity_distribution.bin_edges,"manual",initial_alkalinity_smoothed).normalise();
initial_alkalinity_smoothed_sampler = Geochemistry_Helpers.Sampler(initial_alkalinity_smoothed_distribution,"latin_hypercube");
initial_alkalinity_smoothed_sampler.getSamples(number_of_samples).shuffle();

initial_alkalinity_smoothed_distribution.location = boron_data.age(1);

% initial_alkalinity_smoothed_distribution.plot();

%% We have starting conditions but can refine d11B_sw
% Collate the viable d11B_sw's and make a distribution
initial_d11B_sw = initial_d11B_CO2.collate("boron").collate("d11B_sw").collate("value");
initial_d11B_sw_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:60,initial_d11B_sw).normalise();

% Smooth the distribution to deal with noise - and build a sampler
initial_d11B_sw_smoothed = smooth(initial_d11B_sw_distribution.probabilities,12);
initial_d11B_sw_smoothed_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"manual",initial_d11B_sw_smoothed).normalise();
initial_d11B_sw_smoothed_sampler = Geochemistry_Helpers.Sampler(initial_d11B_sw_smoothed_distribution,"latin_hypercube");
initial_d11B_sw_smoothed_sampler.getSamples(number_of_samples).shuffle();

% We also have a d11B_sw constraint from the range of measured d11B
% So get the minimum and maximum
raw_d11B_minimum = min(boron_data.d11B);
raw_d11B_maximum = max(boron_data.d11B);

% And their uncertainties
raw_d11B_minimum_uncertainty = boron_data.uncertainty(boron_data.d11B==raw_d11B_minimum);
raw_d11B_maximum_uncertainty = boron_data.uncertainty(boron_data.d11B==raw_d11B_maximum);

% Create distributions and samplers
raw_d11B_minimum_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"Gaussian",[raw_d11B_minimum,raw_d11B_minimum_uncertainty]).normalise();
raw_d11B_maximum_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"Gaussian",[raw_d11B_maximum,raw_d11B_maximum_uncertainty]).normalise();
raw_d11B_minimum_sampler = Geochemistry_Helpers.Sampler(raw_d11B_minimum_distribution,"latin_hypercube");
raw_d11B_maximum_sampler = Geochemistry_Helpers.Sampler(raw_d11B_maximum_distribution,"latin_hypercube");

% Get the samples
raw_d11B_minimum_sampler.getSamples(number_of_samples).shuffle();
raw_d11B_maximum_sampler.getSamples(number_of_samples).shuffle();

% The d11B_sw can be anywhere between measured d11B and measured d11B +
% epsilon
% So get random values between 0 and 1 to represent what fraction of
% epsilon is expressed
amount_of_epsilon_distribution = Geochemistry_Helpers.Distribution(-0.1:0.1:1.1,"flat",[0,1]).normalise();
amount_of_epsilon_sampler = Geochemistry_Helpers.Sampler(amount_of_epsilon_distribution,"latin_hypercube");
amount_of_epsilon_sampler.getSamples(number_of_samples).shuffle();

% Calculate estimates of d11B_sw from the minimum
d11B_sw_from_minimum = raw_d11B_minimum_sampler.samples+(amount_of_epsilon_sampler.samples.*epsilon_sampler.samples);

% Shuffle the epsilons so we don't get paired samples
amount_of_epsilon_sampler.shuffle();
% Calculate estimates of d11B_sw from the maximum
d11B_sw_from_maximum = raw_d11B_maximum_sampler.samples+(amount_of_epsilon_sampler.samples.*epsilon_sampler.samples);

% Turn the samples into a distribution
d11B_sw_from_minimum_distribution = Geochemistry_Helpers.Distribution.fromSamples(initial_d11B_sw_distribution.bin_edges,d11B_sw_from_minimum).normalise();
d11B_sw_from_maximum_distribution = Geochemistry_Helpers.Distribution.fromSamples(initial_d11B_sw_distribution.bin_edges,d11B_sw_from_maximum).normalise();

% Multiply them together to get the mutually inclusive region - and put
% that into a distribution
d11B_sw_from_minimum_maximum = d11B_sw_from_minimum_distribution.probabilities.*d11B_sw_from_maximum_distribution.probabilities;
d11B_sw_from_minimum_maximum_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"manual",d11B_sw_from_minimum_maximum).normalise();

% Multiply the d11B_sw from the carbonate chemistry method and d11B
% measured method together for mutually inclusive region
combined_initial_d11B_sw = d11B_sw_from_minimum_maximum_distribution.probabilities.*initial_d11B_sw_smoothed_distribution.probabilities;
combined_initial_d11B_sw_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"manual",combined_initial_d11B_sw).normalise();
combined_initial_d11B_sw_sampler = Geochemistry_Helpers.Sampler(combined_initial_d11B_sw_distribution,"latin_hypercube");
combined_initial_d11B_sw_sampler.getSamples(number_of_samples).shuffle();

% Can do stats on that distribution to see how sure we are about d11B_sw
d11B_sw_mean = mean(combined_initial_d11B_sw_sampler.samples);
d11B_sw_2std = 2*std(combined_initial_d11B_sw_sampler.samples);
d11B_sw_95 = [combined_initial_d11B_sw_distribution.quantile(0.025),combined_initial_d11B_sw_distribution.quantile(0.975)];

% clf
% hold on
% d11B_sw_from_minimum_distribution.plot();
% d11B_sw_from_maximum_distribution.plot();
% d11B_sw_from_minimum_maximum_distribution.plot();
% initial_d11B_sw_smoothed_distribution.plot();
% combined_initial_d11B_sw_distribution.plot();
% 
% xlabel("\delta^{11}B_{sw}");
% ylabel("Probability");
% 
% legend(["From minimum","From maximum","Combined","From carbonate system","Overlap"],'Box','Off','Location','NorthWest');

%% To calculate evolutions we need to know how alkalinity may have changed
% Do this using Gaussian process
gp = Geochemistry_Helpers.GaussianProcess("rbf",boron_data.age');
gp.observations = initial_alkalinity_smoothed_distribution;

gp.runKernel([1700,0.1]);
gp.getSamples(number_of_samples);

% figure(1);
% clf
% hold on;
% for observation_index = 1:numel(gp.observations)
%     plot([gp.observations(observation_index).location,gp.observations(observation_index).location],gp.observations(observation_index).mean()+2*[-gp.observations(observation_index).standard_deviation(),gp.observations(observation_index).standard_deviation()],'Color','k','LineWidth',1);
%     p(3) = plot(gp.observations(observation_index).location,gp.observations(observation_index).mean(),'Color','k','Marker','x');
% end
% plot(gp.queries.collate("location"),gp.queries.quantile(0.5),'Color','k','Marker','.');
% plot(gp.queries.collate("location"),gp.samples(1:10,:));
% 
% clf
% hold on
% h = histogram(gp.samples(:,1),initial_alkalinity_smoothed_distribution.bin_edges,'Normalization','Probability');
% initial_alkalinity_smoothed_distribution.plot();

%% Now that we have alkalinity evolutions, combine with d11B evolution for CO2 evolution
% Create a matrix to hold the data + do the calculations
co2_evolutions = BuCC.d11BCO2().create([height(boron_data),number_of_samples]);
d11B_evolutions = d11B_samplers.collate("samples");
for evolution_index = 1:size(co2_evolutions,2)
    % Assign all the necessary values
    co2_evolutions(:,evolution_index).collate("species_calibration").collate("d11B_measured").assignToEach("value",d11B_evolutions(:,evolution_index));
    
    co2_evolutions(:,evolution_index).collate("boron").assignToAll("epsilon",epsilon_sampler.samples(evolution_index));
    
    count = 0;
    random_index = evolution_index;
    while any(combined_initial_d11B_sw_sampler.samples(random_index)>(d11B_evolutions(:,evolution_index)+epsilon_sampler.samples(evolution_index))) && count<10000;
        random_index = ceil(number_of_samples*rand(1));
        count = count+1;
    end
    co2_evolutions(:,evolution_index).collate("boron").collate("d11B_sw").assignToAll("value",combined_initial_d11B_sw_sampler.samples(random_index));
    
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToEach("temperature",initial_temperature_sampler.samples(evolution_index)+temperature_change);
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("salinity",35);
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("oceanic_pressure",0);
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("atmospheric_pressure",1);
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("calcium",ca_sampler.samples(evolution_index));
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("magnesium",mg_sampler.samples(evolution_index));
    
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").assignToEach("alkalinity",gp.samples(evolution_index,:));
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("equilibrium_coefficients").assignToAll("MyAMI",myami);
end
% Do the calculation
co2_evolutions.calculate();

% Collate the values
d11B_out = co2_evolutions.collate("boron").collate("d11B_4").collate("value");
pH = co2_evolutions.collate("carbonate_chemistry").collate("pH").collate("pValue");
co2 = co2_evolutions.collate("carbonate_chemistry").collate("atmospheric_co2").collate("mole_fraction");
saturation_state = co2_evolutions.collate("carbonate_chemistry").collate("saturation_state");

% There needs to be some postprocessing
% Sometimes (because of the uncertainties on each parameter) you can choose
% a combination that doesn't work mathematically
% That'll give an imaginary result, so remove those
% And if any result in the whole time series was imaginary, then the rest
% aren't valid, so remove them
co2(imag(co2)~=0)=NaN;
for evolution_index = 1:size(co2,2)
    if any(isnan(co2(:,evolution_index)))
        co2(:,evolution_index) = NaN;
    end
end
pH(isnan(co2)) = NaN;
saturation_state(isnan(co2)) = NaN;

%%
pH_change = Geochemistry_Helpers.Distribution.fromSamples((0:0.01:3)',nanmean(pH(1:9,:))-nanmin(pH(10:end,:))).normalise();
% pH_change.plot();

pH_change.quantile(0.5);
pH_change.standard_deviation();

co2_change = Geochemistry_Helpers.Distribution.fromSamples((0:100:10000)',abs(nanmean(co2(1:9,:))-nanmax(co2(10:end,:)))).normalise();
% co2_change.plot();
% xlabel("CO2 (ppm)");
% ylabel("Probability");

%%
% If you want to subsample it can be done like this
pH_subsample_boolean = repmat(saturation_state(1,:)<12 & all(saturation_state<12) & pH(1,:)>8.0,size(pH,1),1);
pH_subsample = reshape(pH(pH_subsample_boolean),22,[]);
co2_subsample = reshape(co2(pH_subsample_boolean),22,[]);
co2_subsample_change =  Geochemistry_Helpers.Distribution.fromSamples((0:100:20000)',abs(nanmean(co2_subsample(1:9,:))-nanmax(co2_subsample(10:end,:)))).normalise();
saturation_state_subsample = reshape(saturation_state(pH_subsample_boolean),22,[]);

co2_subsample_change.quantile(0.5);
co2_subsample_change.standard_deviation()*2;

%
% If you want to subsample it can be done like this
pH_subsample_boolean_2 = repmat(saturation_state(1,:)<12 & all(saturation_state<12)  & pH(1,:)<8.0,size(pH,1),1);
pH_subsample_2 = reshape(pH(pH_subsample_boolean_2),22,[]);
co2_subsample_2 = reshape(co2(pH_subsample_boolean_2),22,[]);
co2_subsample_2_change =  Geochemistry_Helpers.Distribution.fromSamples((0:100:20000)',abs(nanmean(co2_subsample_2(1:9,:))-nanmax(co2_subsample_2(10:end,:))));


co2_subsample_2_change.quantile(0.5);
co2_subsample_2_change.standard_deviation()*2;

%%
% figure(1);
% clf
% hold on
% 
% plot(boron_data.absolute_age,co2_subsample,'b');
% plot(boron_data.absolute_age,co2_subsample_2,'r');
% 
% set(gca,'XDir','Reverse');
% ylim([0,10000]);
% 
% clf
% hold on
% co2_subsample_change.plot();
% co2_subsample_2_change.plot();
% 
% l = legend([">8","<8"]);
% title(l,'Initial pH');
% 
% 
% ax = gca;
% ax.XAxis.Exponent = 0;
% 
% xlabel("\DeltaCO2");
% ylabel("Probability");

toc

%%
% figure(1);
% clf
% hold on
% plot(boron_data.absolute_age,saturation_state_subsample);
% plot(boron_data.absolute_age,median(saturation_state_subsample,2),'k','LineWidth',2);
% set(gca,'XDir','Reverse');

%%

