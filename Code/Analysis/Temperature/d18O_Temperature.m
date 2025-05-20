clear

data_directory = "./../../../Data/";

d18O_sw_datapoints = jsondecode(fileread(data_directory+"/Temperature/d18O_sw_Petryshyn.json")).d18O_sw;
d11B = readtable(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx","Sheet","With_Age");

number_of_samples = 10000;

%% Remove the diagenetic d18O's
d11B_without_diagenetic = d11B;
d11B_without_diagenetic(logical(d11B.("diagenetic_alteration")),:) = [];

% Replace them
d18O_samples = NaN(height(d11B_without_diagenetic),1);
for d18O_index = 1:height(d11B_without_diagenetic)
    d18O_samplers(d18O_index) = Geochemistry_Helpers.Sampler(-10:0.0001:10,"Gaussian",[d11B_without_diagenetic{d18O_index,"d18O"},d11B_without_diagenetic{d18O_index,"d18O_uncertainty"}],'latin_hypercube_random').normalise();
    d18O_samplers(d18O_index).location = d11B_without_diagenetic{d18O_index,"age"};
end

d18O_gp = Geochemistry_Helpers.GaussianProcess("rbf",d11B.("age"));
d18O_gp.observations = d18O_samplers;
d18O_gp.runKernel([1,0.02]);
d18O_gp.getSamples(number_of_samples);

% Fill in estimated d18O's
replaced_d18O_indices = find(d11B.diagenetic_alteration);
for index = replaced_d18O_indices'
    d11B{index,"d18O"} = mean(d18O_gp.samples(:,index));
    d11B{index,"d18O_uncertainty"} = 2*std(d18O_gp.samples(:,index));
end


%% d18Osw
d18O_x = -5:0.01:5;

d18O_sw_samplers = Geochemistry_Helpers.Sampler().create(size(d18O_sw_datapoints,1));
for sampler_index = 1:size(d18O_sw_datapoints,1)
    d18O_sw_samplers(sampler_index) = Geochemistry_Helpers.Sampler(d18O_x,"Gaussian",[d18O_sw_datapoints(sampler_index,1),d18O_sw_datapoints(sampler_index,2)],"latin_hypercube").normalise();
end
d18O_sw_samplers.getSamples(number_of_samples).shuffle();

d18O_sw_combined_sampler = Geochemistry_Helpers.Sampler.fromSamples(d18O_x,mean(d18O_sw_samplers.samples),"latin_hypercube").normalise();
d18O_sw_combined_sampler.getSamples(number_of_samples).shuffle();

d18O_sw = Geochemistry_Helpers.delta().create([number_of_samples,1]);
d18O_sw.assignToAll("standard","O_VSMOW");
d18O_sw.assignToEach("value",d18O_sw_combined_sampler.samples);

%% Anonymous functions
oneil_calibration = @(d18O_ratio,seawater_d18O_ratio) ((((3.597122e-4)*log(d18O_ratio./seawater_d18O_ratio))+1.2194e-6).^(-0.5))-273.15;
kim_oneil_calibration = @(d18O_ratio,seawater_d18O_ratio) ((0.0554631*log(d18O_ratio./seawater_d18O_ratio) + 0.0017981).^(-1))-273.15;
no_ice_hansen_calibration_sw = @(d18O,d18O_sw) -4*(d18O-d18O_sw)+16.8;
anderson_arthur_calibration = @(d18O,d18O_sw) 16-(4.14*(d18O-d18O_sw))+(0.13*(d18O-d18O_sw).^2);

%%
hansen_temperature_matrix = NaN(number_of_samples,height(d11B));
for index = 1:height(d11B)
    d18O_sampler = Geochemistry_Helpers.Sampler(d18O_x,"Gaussian",[d11B.d18O(index),d11B.d18O_uncertainty(index)],"latin_hypercube").normalise();
    d18O_sampler.getSamples(number_of_samples).shuffle();
    
    d18O = Geochemistry_Helpers.delta().create([number_of_samples,1]);
    d18O.assignToAll("standard","O_VPDB");
    d18O.assignToEach("value",d18O_sampler.samples);
    
    oneil_temperatures = oneil_calibration(d18O.collate("ratio"),d18O_sw.collate("ratio"));
    d11B.oneil_temperature(index) = mean(oneil_temperatures);
    d11B.oneil_temperature_uncertainty(index) = 2*std(oneil_temperatures);
    
    kim_oneil_temperatures = kim_oneil_calibration(d18O.collate("ratio"),d18O_sw.collate("ratio"));
    d11B.kim_oneil_temperature(index) = mean(kim_oneil_temperatures);
    d11B.kim_oneil_temperature_uncertainty(index) = 2*std(kim_oneil_temperatures);
    
    hansen_temperatures = no_ice_hansen_calibration_sw(d18O.value,d18O_sw.value);
    d11B.hansen_temperature(index) = mean(hansen_temperatures);
    d11B.hansen_temperature_uncertainty(index) = 2*std(hansen_temperatures);
    
    anderson_arthur_temperatures = anderson_arthur_calibration(d18O.value,d18O_sw.value);
    d11B.anderson_arthur_temperature(index) = mean(anderson_arthur_temperatures);
    d11B.anderson_arthur_temperature_uncertainty(index) = 2*std(anderson_arthur_temperatures);
   
    hansen_temperature_matrix(:,index) = hansen_temperatures;
end
delta_temperatures = hansen_temperature_matrix-hansen_temperature_matrix(:,end);
d11B.delta_temperature = mean(delta_temperatures,1)';
d11B.delta_temperature_uncertainty = 2*std(delta_temperatures,1)';


%% Tiepoint temperature
temperature_x = -50:0.01:50;
combined_samples = [oneil_temperatures,kim_oneil_temperatures,hansen_temperatures,anderson_arthur_temperatures];
combined_sampler = Geochemistry_Helpers.Sampler.fromSamples(temperature_x,combined_samples(:),"latin_hypercube").normalise();

% Fit
temperature_fit_function = @(a1,b1,c1,a2,b2,c2,x) (a1/(c1*sqrt(2*pi)))*exp(-0.5*((x-b1)/c1).^2) + (a2/(c2*sqrt(2*pi)))*exp(-0.5*((x-b2)/c2).^2);
temperature_fit = fit(combined_sampler.bin_midpoints',combined_sampler.probabilities,temperature_fit_function,'StartPoint',[0.05,22,1,0.05,25,1]);

% If you want to check the fit is good:
figure()
hold on
combined_sampler.plot();
plot(temperature_x,temperature_fit(temperature_x));

xlim([20,30]);


%% Save
current_file = data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx";
writematrix(string(d11B.Properties.VariableNames),current_file,"Sheet","Temperature_Calibrations","Range","A1");
writematrix([" ","m","Ma","‰","‰","‰","‰","‰","‰","","","°C","°C","°C","°C","°C"],current_file,"Sheet","Temperature_Calibrations","Range","A2");
writecell(d11B{:,1},current_file,"Sheet","Temperature_Calibrations","Range","A3");
writematrix(d11B{:,2:end},current_file,"Sheet","Temperature_Calibrations","Range","B3");

% Tiepoint
pH_change_parameters = jsondecode(fileread(data_directory+"/pH_Change/Input.json"));
pH_change_parameters.tiepoint_temperature = round([temperature_fit.a1,temperature_fit.b1,temperature_fit.c1,temperature_fit.a2,temperature_fit.b2,temperature_fit.c2],3);

writeInputJSON(data_directory+"/pH_Change/Input.json",pH_change_parameters);

