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
% Choose the number of statistical samples
number_of_samples = 1000;

% Load in the data
temperature_data = readtable("./../Data/TJ_Temperature.xlsx");
boron_data = readtable("./../Data/TJ_d11B_pH.xlsx");

boron_data.age = boron_data.Absolute_age_Ma;

% Segregate into before + during/after perturbation
background_data = boron_data(1:9,:);
background_d11B4 = nanmean(background_data.d11B);
background_d11B4_with_uncertainty = background_d11B4+normrnd(0,0.5,number_of_samples,1);

perturbation_data = boron_data(10:end,:);
perturbation_d11B4 = nanmin(perturbation_data.d11B);

% Generate distributions for each variable to represent uncertainties
initial_co2_distribution = Geochemistry_Helpers.Distribution(0:100e-6:10000e-6,"Flat",[500e-6,5000e-6]).normalise();
initial_omega_distribution = Geochemistry_Helpers.Distribution(0:0.1:12,"Flat",[5,10.7]).normalise();
initial_temperature_distribution = Geochemistry_Helpers.Distribution(-10:5:60,"Flat",[10,30]).normalise();

% Use a sampling technique for those distributions
initial_co2_sampler = Geochemistry_Helpers.Sampler(initial_co2_distribution,"latin_hypercube");
initial_omega_sampler = Geochemistry_Helpers.Sampler(initial_omega_distribution,"latin_hypercube");
initial_temperature_sampler = Geochemistry_Helpers.Sampler(initial_temperature_distribution,"latin_hypercube");

% Get the samples
initial_co2_sampler.getSamples(number_of_samples).shuffle();
initial_omega_sampler.getSamples(number_of_samples).shuffle();
initial_temperature_sampler.getSamples(number_of_samples).shuffle();

% Create an array for carbonate chemistry calculations
initial_carbonate_chemistry = BuCC.CarbonateChemistry().create(number_of_samples);

% Assign the parameters
initial_carbonate_chemistry.collate("atmospheric_co2").assignToEach("partial_pressure",initial_co2_sampler.samples);
initial_carbonate_chemistry.assignToEach("saturation_state",initial_omega_sampler.samples);

initial_carbonate_chemistry.assignToEach("temperature",initial_temperature_sampler.samples);

initial_carbonate_chemistry.assignToAll("units"," mol/kg");
initial_carbonate_chemistry.assignToAll("salinity",35);
initial_carbonate_chemistry.assignToAll("oceanic_pressure",0);
initial_carbonate_chemistry.assignToAll("atmospheric_pressure",1);
initial_carbonate_chemistry.assignToAll("calcium",17);
initial_carbonate_chemistry.assignToAll("magnesium",20);

% Create a MyAMI object
myami = MyAMI.MyAMI("Precalculated",true);
initial_carbonate_chemistry.collate("equilibrium_coefficients").assignToAll("MyAMI",myami);

% Do carbonate chemistry calculations
initial_carbonate_chemistry.calculate();

% Collate the results of those calculations back to normal arrays
initial_pH = initial_carbonate_chemistry.collate("pH").collate("pValue");
initial_alkalinity = initial_carbonate_chemistry.collate("alkalinity")*1e6; %(converted units from mol/kg to umol/kg)

% Create a distribution from those sample arrays
initial_pH_distribution = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,initial_pH).normalise();
initial_alkalinity_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:100:12000,initial_alkalinity).normalise();

% Smooth the distribution of alkalinity to minimise noise (from statistical
% randomness)
initial_alkalinity_smoothed = smooth(initial_alkalinity_distribution.probabilities,12);
% Create a new distribution from the smoothed data and a sampler to get
% plausible alkalinities
initial_alkalinity_smoothed_distribution = Geochemistry_Helpers.Distribution(initial_alkalinity_distribution.bin_edges,"manual",initial_alkalinity_smoothed).normalise();
initial_alkalinity_smoothed_sampler = Geochemistry_Helpers.Sampler(initial_alkalinity_smoothed_distribution,"latin_hypercube");
initial_alkalinity_smoothed_sampler.getSamples(number_of_samples).shuffle();

% clf
% hold on
% initial_alkalinity_distribution.plot();
% initial_alkalinity_smoothed_distribution.plot();
% 
% xlabel("Alkalinity \mumol/kg");
% ylabel("Probability");

%% Now we have pH do the boron calculations
% Create uncertainty on epsilon
epsilon_possibilities = Geochemistry_Helpers.Distribution([24:0.01:31],"Gaussian",[27.2,0.6]).normalise();
epsilon_possibilities_sampler = Geochemistry_Helpers.Sampler(epsilon_possibilities,"latin_hypercube");
epsilon_possibilities_sampler.getSamples(number_of_samples).shuffle();

% Create an array to do boron calculations
initial_boron = BuCC.BoronpH().create(number_of_samples);
% Assign the data (pH, measured d11B, epsilon, T/S/P/Ca/Mg for Kb)
initial_boron.collate("pH").assignToEach("pValue",initial_pH);
initial_boron.collate("d11B_4").assignToEach("value",background_d11B4_with_uncertainty);
initial_boron.assignToEach("epsilon",epsilon_possibilities_sampler.samples);

conditions = initial_boron.collate("pKb").collate("conditions");
conditions.assignToAll("temperature",25);
conditions.assignToAll("salinity",35);
conditions.assignToAll("oceanic_pressure",0);
conditions.assignToAll("atmospheric_pressure",1);
conditions.assignToAll("calcium",17);
conditions.assignToAll("magnesium",20);

% Calculate pKb
initial_boron.collate("pKb").calculate();

% Make sure d11B_sw is clear
initial_boron.collate("d11B_sw").assignToAll("value",NaN);

% Do the boron calculations
initial_boron.calculate();

%% We have starting conditions but can refine d11B_sw
% Collate the viable d11B_sw's and make a distribution
initial_d11B_sw = initial_boron.collate("d11B_sw").collate("value");
initial_d11B_sw_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:40,initial_d11B_sw).normalise();

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
raw_d11B_minimum_uncertainty = boron_data.x2SE(boron_data.d11B==raw_d11B_minimum);
raw_d11B_maximum_uncertainty = boron_data.x2SE(boron_data.d11B==raw_d11B_maximum);

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
epsilon_distribution = Geochemistry_Helpers.Distribution(-0.1:0.1:1.1,"flat",[0,1]).normalise();
epsilon_sampler = Geochemistry_Helpers.Sampler(epsilon_distribution,"latin_hypercube");
epsilon_sampler.getSamples(number_of_samples).shuffle();

% Calculate estimates of d11B_sw from the minimum
d11B_sw_from_minimum = raw_d11B_minimum_sampler.samples+(epsilon_sampler.samples*27.2);

% Shuffle the epsilons so we don't get paired samples
epsilon_sampler.shuffle();
% Calculate estimates of d11B_sw from the maximum
d11B_sw_from_maximum = raw_d11B_maximum_sampler.samples+(epsilon_sampler.samples*27.2);

% Turn the samples into a distribution
d11B_sw_from_minimum_distribution = Geochemistry_Helpers.Distribution.fromSamples(initial_d11B_sw_distribution.bin_edges,d11B_sw_from_minimum).normalise();
d11B_sw_from_maximum_distribution = Geochemistry_Helpers.Distribution.fromSamples(initial_d11B_sw_distribution.bin_edges,d11B_sw_from_maximum).normalise();

% Multiply them together to get the mutually inclusive region - and put
% that into a distribution
d11B_sw_from_minimum_maximum = d11B_sw_from_minimum_distribution.probabilities.*d11B_sw_from_maximum_distribution.probabilities;
d11B_sw_from_minimum_maximum_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"manual",d11B_sw_from_minimum_maximum).normalise();

% Multiply the d11B_sw from the carbonate chemistry method and d11B
% measured method together for mutually inclusive region
combined_initial_d11B_sw = d11B_sw_from_minimum_maximum_distribution.probabilities.*initial_d11B_sw_smoothed_distribution.probabilities';
combined_initial_d11B_sw_distribution = Geochemistry_Helpers.Distribution(initial_d11B_sw_distribution.bin_edges,"manual",combined_initial_d11B_sw).normalise();
combined_initial_d11B_sw_sampler = Geochemistry_Helpers.Sampler(combined_initial_d11B_sw_distribution,"latin_hypercube");
combined_initial_d11B_sw_sampler.getSamples(number_of_samples).shuffle();

% Can do stats on that distribution to see how sure we are about d11B_sw
d11B_sw_mean = mean(combined_initial_d11B_sw_sampler.samples);
d11B_sw_2std = 2*std(combined_initial_d11B_sw_sampler.samples);
d11B_sw_95 = [combined_initial_d11B_sw_distribution.quantile(0.025),combined_initial_d11B_sw_distribution.quantile(0.975)];

%
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
% Put the mean alkalinity into the table
boron_data.alkalinity = [mean(initial_alkalinity_smoothed_sampler.samples);NaN(height(boron_data)-1,1)];
% We only have constraints on the starting value
boron_observations = boron_data(1,:);

age_midpoints = boron_data.age';
age_matrix = repmat(age_midpoints,number_of_samples,1);

% Define an anonymous function for kernel - radial basis function
kernel_function = @(t,t1,s,l) (s^2)*exp(-((norm(t-t1)/l).^2)/2);

% Two parameters to control variation in each axis
vertical_scaling = 1200; %(in umol/kg of alkalinity)
horizontal_scaling = 0.1; %(in Myr)

% Run the kernel for the age matrix - tells you how related two times are
k_oo = runKernel(kernel_function,boron_observations.age,boron_observations.age,vertical_scaling,horizontal_scaling)';
k_os = runKernel(kernel_function,boron_observations.age,age_matrix(1,:),vertical_scaling,horizontal_scaling)';
k_so = k_os'; %runKernel(kernel_function,age_matrix(1,:),boron_observations.age,vertical_scaling,horizontal_scaling)';
k_ss = runKernel(kernel_function,age_matrix(1,:),age_matrix(1,:),vertical_scaling,horizontal_scaling)';

k_oo = k_oo+diag(repelem(1e-6,size(k_oo,1))); % Add a very small amount for numerical stability

% We started with a single value for alkalinity so need to add uncertainty
% to that
observation_uncertainty = ((2*var(initial_alkalinity_smoothed_sampler.samples))).*eye(size(k_oo));

% Matrix maths for calculation the mean at each place and uncertainty
interpolation_mean = nanmean(boron_observations.alkalinity) + k_so/(k_oo+observation_uncertainty)*(boron_observations.alkalinity-nanmean(boron_observations.alkalinity));
interpolation_sigma = k_ss-k_so/(k_oo+observation_uncertainty)*k_os;

% Ensure the matrix is upper triangular (this is just for numerical
% stability)
sigma_difference = triu(interpolation_sigma-interpolation_sigma');
interpolation_sigma = interpolation_sigma-sigma_difference;
interpolation_sigma = interpolation_sigma + max(-eig(interpolation_sigma))*speye(size(interpolation_sigma));

% Generate random samples using the calculated parameters
interpolated_alkalinity = mvnrnd(interpolation_mean,interpolation_sigma,1000);

% clf
% plot(age_matrix(1,:),interpolated_alkalinity(1:1000,:));
% set(gca,'XDir','Reverse');
% ylim([1000,9000]);
% 
% xlabel("Age (Ma)");
% ylabel("Alkalinity (\mumol/kg)");


%% Now that we have alkalinity evolutions, combine with d11B evolution for CO2 evolution
% Create a matrix to hold the data + do the calculations
co2_evolutions = BuCC.d11BCO2().create([height(boron_data),1000]);
for evolution_index = 1:size(co2_evolutions,2)
    % Assign all the necessary values
    co2_evolutions(:,evolution_index).collate("species_calibration").collate("d11B_measured").assignToEach("value",boron_data.d11B);
    co2_evolutions(:,evolution_index).collate("boron").collate("d11B_sw").assignToAll("value",combined_initial_d11B_sw_sampler.samples(evolution_index));
    co2_evolutions(:,evolution_index).collate("boron").assignToAll("epsilon",27.2);
    
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("temperature",25);
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("salinity",35);
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("oceanic_pressure",0);
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("atmospheric_pressure",0);
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("calcium",17);
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("conditions").assignToAll("magnesium",20);
    
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").assignToEach("alkalinity",interpolated_alkalinity(evolution_index,:));
    co2_evolutions(:,evolution_index).collate("carbonate_chemistry").collate("equilibrium_coefficients").assignToAll("MyAMI",myami);
end
% Do the calculation
co2_evolutions.calculate();

% Collate the values
pH = co2_evolutions.collate("carbonate_chemistry").collate("pH").collate("pValue");
co2 = co2_evolutions.collate("carbonate_chemistry").collate("atmospheric_co2").collate("partial_pressure");

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

% If you want to subsample it can be done like this
pH_subsample_boolean = repmat(pH(1,:)>7.95 & pH(1,:)<8,size(pH,1),1);
pH_subsample = reshape(pH(pH_subsample_boolean),22,[]);
co2_subsample = reshape(co2(pH_subsample_boolean),22,[]);

%%
function covariance_matrix = runKernel(kernel_function,x_1,x_2,sigma,l)
    for x_1_index = 1:numel(x_1)
        for x_2_index = 1:numel(x_2)
            covariance_matrix(x_2_index,x_1_index) = kernel_function(x_1(x_1_index),x_2(x_2_index),sigma,l);
        end
    end
end