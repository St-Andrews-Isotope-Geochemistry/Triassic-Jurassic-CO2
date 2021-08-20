clear

data_directory = "./../../Data/";
d18O_d13C = readtable(data_directory+"TJ_d18O_d13C.xlsx","Sheet","With_Age");

%%

% anderson_arthur_calibration = @(d18O,d18O_sw) 16-(4.14*(d18O-d18O_sw))+(0.13*(d18O-d18O_sw).^2);
% no_ice_hansen_calibration = @(d18O) -4*d18O+c; @ d18O_sw = -1.2;
% If it was originally d18O-d18O_sw as in arthur-anderson:
% -4*(d18O-d18O_sw)+12 -> -4*+1.2 = -4.8 -> -4*(d18O-d18O_sw)+16.8

no_ice_hansen_calibration_sw = @(d18O,d18O_sw) -4*(d18O-d18O_sw)+16.8;

%%
d18O_sw = [2.1,2.3,2.6,-0.6,2.7,2.4,3.0,3.3,2.3,-1.3,-1.3,0,-1.6,-1.2];
d18O_sw_uncertainty = [1.1,0.9,0.9,1.1,1.0,0.9,0.8,0.7,0.7,0.8,1.0,0.7,0.7,0.6];

d18O_x = -5:0.01:5;
number_of_samples = 1000;

for sampler_index = 1:numel(d18O_sw)
    d18O_sw_samplers(sampler_index) = Geochemistry_Helpers.Sampler(d18O_x,"Gaussian",[d18O_sw(sampler_index),d18O_sw_uncertainty(sampler_index)],"latin_hypercube").normalise();
end
d18O_sw_samplers.getSamples(number_of_samples).shuffle();

d18O_sw_combined_distribution = Geochemistry_Helpers.Distribution.fromSamples(d18O_x,mean(d18O_sw_samplers.samples));
d18O_sw_95 = d18O_sw_combined_distribution.quantile([0.025,0.5,0.975]);
d18O_sw_metrics = [d18O_sw_combined_distribution.median(),d18O_sw_combined_distribution.standard_deviation()];

%%
d18O_d13C.hansen_temperature = no_ice_hansen_calibration_sw(d18O_d13C.d18O,d18O_sw_metrics(1));

%%
clf
plot(d18O_d13C.age,d18O_d13C.hansen_temperature,'.');
set(gca,'XDir','Reverse');