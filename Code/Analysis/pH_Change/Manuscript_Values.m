clear

%%
data_directory = "./../../../Data/";

%% Background d11B + Minimum Change
% age = jsondecode(fileread(data_directory+"/pH_change/prior.json")).age;
prior = jsondecode(fileread(data_directory+"/pH_change/prior.json")).prior;
posterior = jsondecode(fileread(data_directory+"/pH_change/posterior.json")).posterior;

%% Pre excursion d11B

%% pH Change
posterior_pH = [posterior.pH];
posterior_initial_pH = posterior_pH(end,:);
posterior_minimum_pH = min(posterior_pH);

posterior_pH_change = -(posterior_initial_pH-posterior_minimum_pH);
posterior_pH_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(-2:0.01:2,posterior_pH_change).normalise();

posterior_pH_change_0975 = posterior_pH_change_distribution.quantile(0.975);
posterior_pH_change_0950 = posterior_pH_change_distribution.quantile(0.950);
posterior_pH_change_0500 = posterior_pH_change_distribution.quantile(0.50);
posterior_pH_change_0050 = posterior_pH_change_distribution.quantile(0.050);
posterior_pH_change_0025 = posterior_pH_change_distribution.quantile(0.025);

%% CO2 Change
posterior_co2 = [posterior.co2];
posterior_initial_co2 = posterior_co2(end,:);
posterior_maximum_co2 = max(posterior_co2);

posterior_co2_change = posterior_maximum_co2-posterior_initial_co2;
posterior_co2_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(-2000:1000:100000,posterior_co2_change).normalise();

posterior_co2_change_005 = posterior_co2_change_distribution.quantile(0.05);

posterior_co2_doublings = log2(posterior_maximum_co2)-log2(posterior_initial_co2);
posterior_co2_doublings_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:10,posterior_co2_doublings).normalise();
posterior_co2_doublings_quantiles = posterior_co2_doublings_distribution.quantile([0.025,0.05,0.05,0.95,0.975]);

%% Initial CO2
posterior_initial_co2 = [posterior.co2_initial];
posterior_co2_initial_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:100:100000,posterior_initial_co2).normalise();

posterior_co2_initial_quantiles = posterior_co2_initial_distribution.quantile([0.025,0.975]);

posterior_co2 = [posterior.co2];
delta_co2 = log2(posterior_co2)-log2(posterior_co2(end,:));


%% Species calibration gradient
posterior_species_calibration = [posterior.species_calibration];
posterior_species_calibration_gradient = posterior_species_calibration(1,:);
posterior_species_calibration_gradient_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:1,posterior_species_calibration_gradient).normalise();

posterior_species_calibration_gradient_005 = posterior_species_calibration_gradient_distribution.quantile(0.05);

%%
if exist(data_directory+"/pH_Change/Manuscript_Values_standard_smooth.json",'file')
    existing_manuscript_values = jsondecode(fileread(data_directory+"/pH_Change/Manuscript_Values_standard_smooth.json"));
end
existing_manuscript_values.pH_change_0025_quantile = round(posterior_pH_change_0975,3);
existing_manuscript_values.pH_change_0050_quantile = round(posterior_pH_change_0950,3);
existing_manuscript_values.pH_change_0500_quantile = round(posterior_pH_change_0500,3);
existing_manuscript_values.pH_change_0950_quantile = round(posterior_pH_change_0050,3);
existing_manuscript_values.pH_change_0975_quantile = round(posterior_pH_change_0025,3);
existing_manuscript_values.co2_change_005_quantile = round(posterior_co2_change_005,3);
existing_manuscript_values.initial_co2 = round(posterior_co2_initial_quantiles,3);
existing_manuscript_values.species_calibration_gradient = round(posterior_species_calibration_gradient_005,3);
existing_manuscript_values.co2_doublings_0025_quantile = round(posterior_co2_doublings_quantiles(1),3);
existing_manuscript_values.co2_doublings_0050_quantile = round(posterior_co2_doublings_quantiles(2),3);
existing_manuscript_values.co2_doublings_0500_quantile = round(posterior_co2_doublings_quantiles(3),3);
existing_manuscript_values.co2_doublings_0950_quantile = round(posterior_co2_doublings_quantiles(4),3);
existing_manuscript_values.co2_doublings_0975_quantile = round(posterior_co2_doublings_quantiles(5),3);

writeInputJSON(data_directory+"/pH_Change/Manuscript_Values_standard_smooth.json",existing_manuscript_values);
