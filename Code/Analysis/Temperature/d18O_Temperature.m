clear

data_directory = "./../../../Data/";
% d18O_d13C = readtable(data_directory+"/Temperature/TJ_d18O_d13C.xlsx","Sheet","With_Age");
d18O_sw = jsondecode(fileread(data_directory+"/Temperature/d18O_sw_Petryshyn.json")).d18O_sw;

d11B = readtable(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx","Sheet","With_Age");

%%
d18O_x = -5:0.01:5;
number_of_samples = 10000;

d18O_sw_samplers = Geochemistry_Helpers.Sampler().create(size(d18O_sw,1));
for sampler_index = 1:size(d18O_sw,1)
    d18O_sw_samplers(sampler_index) = Geochemistry_Helpers.Sampler(d18O_x,"Gaussian",[d18O_sw(sampler_index,1),d18O_sw(sampler_index,2)],"latin_hypercube").normalise();
end
d18O_sw_samplers.getSamples(number_of_samples).shuffle();

d18O_sw_combined_distribution = Geochemistry_Helpers.Distribution.fromSamples(d18O_x,mean(d18O_sw_samplers.samples));
d18O_sw_metrics = [d18O_sw_combined_distribution.median(),d18O_sw_combined_distribution.standard_deviation()];

%% Anonymous functions
oneil_calibration = @(d18O_ratio,seawater_d18O_ratio) ((((3.597122e-4)*log(d18O_ratio./seawater_d18O_ratio))+1.2194e-6).^(-0.5))-273.15;
kim_oneil_calibration = @(d18O_ratio,seawater_d18O_ratio) ((0.0554631*log(d18O_ratio./seawater_d18O_ratio) + 0.0017981).^(-1))-273.15;
no_ice_hansen_calibration_sw = @(d18O,d18O_sw) -4*(d18O-d18O_sw)+16.8;
anderson_arthur_calibration = @(d18O,d18O_sw) 16-(4.14*(d18O-d18O_sw))+(0.13*(d18O-d18O_sw).^2);

%%
d18O = Geochemistry_Helpers.delta().create([height(d11B),1]);
d18O.assignToAll("standard","O_VPDB");
d18O.assignToEach("value",d11B.d18O);

d18O_lower = Geochemistry_Helpers.delta().create([height(d11B),1]);
d18O_lower.assignToAll("standard","O_VPDB");
d18O_lower.assignToEach("value",d11B.d18O-d11B.d18O_uncertainty);

d18O_upper = Geochemistry_Helpers.delta().create([height(d11B),1]);
d18O_upper.assignToAll("standard","O_VPDB");
d18O_upper.assignToEach("value",d11B.d18O+d11B.d18O_uncertainty);

seawater_d18O_Petryshyn = Geochemistry_Helpers.delta("O_VSMOW",d18O_sw_metrics(1));

d11B.oneil_temperature = oneil_calibration(d18O.collate("ratio"),seawater_d18O_Petryshyn.ratio);
d11B.kim_oneil_temperature = kim_oneil_calibration(d18O.collate("ratio"),seawater_d18O_Petryshyn.ratio);
d11B.hansen_temperature = no_ice_hansen_calibration_sw(d11B.d18O,d18O_sw_metrics(1));
d11B.anderson_arthur_temperature = anderson_arthur_calibration(d11B.d18O,seawater_d18O_Petryshyn.value);

d11B.oneil_temperature_uncertainty = abs(oneil_calibration(d18O_upper.collate("ratio"),seawater_d18O_Petryshyn.ratio)-oneil_calibration(d18O_lower.collate("ratio"),seawater_d18O_Petryshyn.ratio));
d11B.kim_oneil_temperature_uncertainty = abs(kim_oneil_calibration(d18O_upper.collate("ratio"),seawater_d18O_Petryshyn.ratio)-kim_oneil_calibration(d18O_lower.collate("ratio"),seawater_d18O_Petryshyn.ratio));
d11B.hansen_temperature_uncertainty = abs(no_ice_hansen_calibration_sw(d11B.d18O+d11B.d18O_uncertainty,d18O_sw_metrics(1))-no_ice_hansen_calibration_sw(d11B.d18O-d11B.d18O_uncertainty,d18O_sw_metrics(1)));
d11B.anderson_arthur_temperature_uncertainty = abs(anderson_arthur_calibration(d11B.d18O-d11B.d18O_uncertainty,seawater_d18O_Petryshyn.value)-anderson_arthur_calibration(d11B.d18O+d11B.d18O_uncertainty,seawater_d18O_Petryshyn.value));

background_data = d11B(d11B.age>=d11B.age(9),:);

d11B.delta_temperature = d11B.hansen_temperature-mean(background_data.hansen_temperature,1,"omitnan");


%% Save
current_file = data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx";
writematrix(string(d11B.Properties.VariableNames),current_file,"Sheet","Temperature_Calibrations","Range","A1");
writematrix([" ","m","Ma","‰","‰","‰","‰","‰","‰","°C","°C","°C","°C","°C"],current_file,"Sheet","Temperature_Calibrations","Range","A2");
writecell(d11B{:,1},current_file,"Sheet","Temperature_Calibrations","Range","A3");
writematrix(d11B{:,2:end},current_file,"Sheet","Temperature_Calibrations","Range","B3");

%%
% averaged.height = unique(d18O_d13C.height);
% for depth_index = 1:numel(averaged.height)
%     at_depth = d18O_d13C.height==averaged.height(depth_index);
%     
%     averaged.age(depth_index,1) = mean(d18O_d13C.age(at_depth),'omitnan');
%     
%     averaged.d13C(depth_index,1) = mean(d18O_d13C.d13C(at_depth),'omitnan');
%     averaged.d13C_uncertainty(depth_index,1) = std(d18O_d13C.d13C(at_depth),'omitnan');
%     
%     averaged.d18O(depth_index,1) = mean(d18O_d13C.d18O(at_depth),'omitnan');
%     averaged.d18O_uncertainty(depth_index,1) = std(d18O_d13C.d18O(at_depth),'omitnan');
%     
%     averaged.oneil_temperature(depth_index,1) = mean(d18O_d13C.oneil_temperature(at_depth),'omitnan');
%     averaged.oneil_temperature_uncertainty(depth_index,1) = std(d18O_d13C.oneil_temperature(at_depth),'omitnan')./sqrt(numel(d18O_d13C.oneil_temperature(at_depth)));
%     
%     averaged.kim_oneil_temperature(depth_index,1) = mean(d18O_d13C.kim_oneil_temperature(at_depth),'omitnan');
%     averaged.kim_oneil_temperature_uncertainty(depth_index,1) = std(d18O_d13C.kim_oneil_temperature(at_depth),'omitnan')./sqrt(numel(d18O_d13C.kim_oneil_temperature(at_depth)));
%     
%     averaged.hansen_temperature(depth_index,1) = mean(d18O_d13C.hansen_temperature(at_depth),'omitnan');
%     averaged.hansen_temperature_uncertainty(depth_index,1) = std(d18O_d13C.hansen_temperature(at_depth),'omitnan')./sqrt(numel(d18O_d13C.hansen_temperature(at_depth)));
% 
%     averaged.anderson_arthur_temperature(depth_index,1) = mean(d18O_d13C.anderson_arthur_temperature(at_depth),'omitnan');
%     averaged.anderson_arthur_temperature_uncertainty(depth_index,1) = std(d18O_d13C.anderson_arthur_temperature(at_depth),'omitnan')./sqrt(numel(d18O_d13C.anderson_arthur_temperature(at_depth)));
% end
% d18O_uncertainty_assumed = max(averaged.d18O_uncertainty(averaged.d18O_uncertainty~=0),[],'omitnan');
% d13C_uncertainty_assumed = max(averaged.d13C_uncertainty(averaged.d13C_uncertainty~=0),[],'omitnan');
% oneil_temperature_uncertainty_assumed = max(averaged.oneil_temperature_uncertainty(averaged.oneil_temperature_uncertainty~=0),[],'omitnan');
% kim_oneil_temperature_uncertainty_assumed = max(averaged.kim_oneil_temperature_uncertainty(averaged.kim_oneil_temperature_uncertainty~=0),[],'omitnan');
% hansen_temperature_uncertainty_assumed = max(averaged.hansen_temperature_uncertainty(averaged.hansen_temperature_uncertainty~=0),[],'omitnan');
% anderson_arthur_temperature_uncertainty_assumed = max(averaged.anderson_arthur_temperature_uncertainty(averaged.anderson_arthur_temperature_uncertainty~=0),[],'omitnan');
% 
% averaged.d18O_uncertainty(averaged.d18O_uncertainty==0) = d18O_uncertainty_assumed;
% averaged.d13C_uncertainty(averaged.d13C_uncertainty==0) = d13C_uncertainty_assumed;
% averaged.oneil_temperature_uncertainty(averaged.oneil_temperature_uncertainty==0) = oneil_temperature_uncertainty_assumed;
% averaged.kim_oneil_temperature_uncertainty(averaged.kim_oneil_temperature_uncertainty==0) = kim_oneil_temperature_uncertainty_assumed;
% averaged.hansen_temperature_uncertainty(averaged.hansen_temperature_uncertainty==0) = hansen_temperature_uncertainty_assumed;
% averaged.anderson_arthur_temperature_uncertainty(averaged.anderson_arthur_temperature_uncertainty==0) = anderson_arthur_temperature_uncertainty_assumed;
% 
% averaged = struct2table(averaged);

%% Save
% current_file = data_directory+"/Temperature/TJ_d18O_d13C.xlsx";
% writematrix(string(averaged.Properties.VariableNames),current_file,"Sheet","Averaged","Range","A1");
% writematrix(["m","Ma","‰","‰","‰","‰","°C","°C","°C","°C","°C","°C","°C","°C"],current_file,"Sheet","Averaged","Range","A2");
% writematrix(averaged{:,:},current_file,"Sheet","Averaged","Range","A3");

%% Delta temperature
% delta_temperature = d11B.hansen_temperature-d11B.hansen_temperature(1);
% 
% current_file = data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx";
% writematrix(["height","age","d13C","d13C_uncertainty","d18O","d18O_uncertainty","delta_temperature","delta_temperature_uncertainty"],current_file,"Sheet","Delta_Temperature","Range","A1");
% writematrix(["m","Ma","‰","‰","‰","‰","°C","°C"],current_file,"Sheet","Delta_Temperature","Range","A2");
% writematrix([averaged.height,averaged.age,averaged.d13C,averaged.d13C_uncertainty,averaged.d18O,averaged.d18O_uncertainty,delta_temperature,averaged.hansen_temperature_uncertainty],current_file,"Sheet","Delta_Temperature","Range","A3");

