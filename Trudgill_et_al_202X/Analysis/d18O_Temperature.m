clear
data = readtable("./../Data/d18O_d13C","Sheet","Matlab");

%% Anonymous functions
no_ice_hansen_calibration = @(d18O) -4*d18O+12;
oneil_calibration = @(d18O_ratio,seawater_d18O_ratio) ((((3.597122e-4)*log(d18O_ratio./seawater_d18O_ratio))+1.2194e-6).^(-0.5))-273.15;

%%
d18O(height(data),1) = Geochemistry_Helpers.delta();
d18O.assignToAll("standard","O");

% d18O is measured on the VPDB standard
% But for O'Neil it must be on VSMOW standard
% VPDB = 0.97001×VSMOW−29.99
% So VSMOW = (VPDB+29.99)/0.97001
% Or VSMOW = 1.0309(VPDB+29.99)
% Or VSMOW = 1.0309*VPDB + 30.9172

d18O.assignToEach("value",(1.0309*data.d18O)+30.9172);
seawater_d18O = delta("O",-1.2);

% In O'Neil 1969
% For CaCO3-H2O
% 1000*ln(a) = 2.78(10^6 * T^(-2))-3.39
% Where a = (18O/16O CaCO3)/(18O/16O H2O)
% T = Temperature (in degK)

% So for temperature:
% ((1000*ln(a))+3.39)/2.78 = 10^6 * T^(-2)
% (((1000*ln(a))+3.39)/2.78) / 10^6 = T^(-2)
% ((((1000*ln(a))+3.39)/2.78) / 10^6)^(-0.5) = T;
% Alternatively: ((3.597122e-4 *log(a))+1.2194e-6).^(-0.5)'XDir','Reverse');

% Test case
% Paper says T=25degC at 1000*ln(a)=27.9
% 27.9 = 2.78(10^6 * T^(-2))-3.39
% (((27.9+3.39)/2.78)/(10^6)).^(-0.5) = 298.0709
% Therefore assume in Kelvin not Celcius
% (((27.9+3.39)/2.78)/(10^6)).^(-0.5) - 273.15 = 24.9209
% Mismatch of ~0.08 degC in this case
% To double check my case
% ((((3.597122e-4)*(27.9/1000))+1.2194e-6).^(-0.5))-273.15 = 24.9213
% Difference is likely due to numerical rounding (and this version is
% slightly closer to the true answer).

data.oneil_temperature = oneil_calibration(d18O.collate("ratio"),seawater_d18O.ratio);
data.hansen_temperature = no_ice_hansen_calibration(data.d18O);

%%
averaged.depths = unique(data.Depth_Catherine);
for depth_index = 1:numel(averaged.depths)
    at_depth = data.Depth_Catherine==averaged.depths(depth_index);
    averaged.oneil_temperature(depth_index,1) = nanmean(data.oneil_temperature(at_depth));
    averaged.oneil_temperature_uncertainty(depth_index,1) = nanstd(data.oneil_temperature(at_depth));
    
    averaged.hansen_temperature(depth_index,1) = nanmean(data.hansen_temperature(at_depth));
    averaged.hansen_temperature_uncertainty(depth_index,1) = nanstd(data.hansen_temperature(at_depth));
    
    averaged.age(depth_index,1) = nanmean(data.Absolute_Age_Ma(at_depth));
    averaged.d13C(depth_index,1) = nanmean(data.d13C(at_depth));
    averaged.d18O(depth_index,1) = nanmean(data.d18O(at_depth));
end
averaged.oneil_temperature_uncertainty(averaged.oneil_temperature_uncertainty==0) = nanmean(averaged.oneil_temperature_uncertainty(averaged.oneil_temperature_uncertainty~=0));
averaged.hansen_temperature_uncertainty(averaged.hansen_temperature_uncertainty==0) = nanmean(averaged.hansen_temperature_uncertainty(averaged.hansen_temperature_uncertainty~=0));

averaged_table = struct2table(averaged);

%%
parameters_to_save = ["depths","age","d13C","d18O","hansen_temperature","hansen_temperature_uncertainty"];
parameters_to_save_as = ["depth","age","d13C","d18O","temperature","temperature_uncertainty"];

for parameter_index = 1:numel(parameters_to_save)
    TJ_temperature.(parameters_to_save_as(parameter_index)) = averaged_table.(parameters_to_save(parameter_index));
end

writetable(struct2table(TJ_temperature),"./Data/TJ_Temperature.xlsx");

%%
clf
hold on
for depth_index = 1:height(averaged_table)
    plot([averaged.age(depth_index),averaged.age(depth_index)],averaged.oneil_temperature(depth_index)+[-2*averaged.oneil_temperature_uncertainty(depth_index),2*averaged.oneil_temperature_uncertainty(depth_index)],'r');
    plot([averaged.age(depth_index),averaged.age(depth_index)],averaged.hansen_temperature(depth_index)+[-2*averaged.hansen_temperature_uncertainty(depth_index),2*averaged.hansen_temperature_uncertainty(depth_index)],'b');
end
oneil_handle = plot(averaged.age,averaged.oneil_temperature,'xr');
hansen_handle = plot(averaged.age,averaged.hansen_temperature,'xb');
xlabel("Age");
ylabel("Temperature ^{\circ}C");

set(gca,"XDir","Reverse");

legend_handle = legend([oneil_handle,hansen_handle],["O-Neil","Hansen"],"Location","NorthWest","Box","Off");