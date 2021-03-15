clear
data = readtable("./../Data/d18O_d13C","Sheet","Matlab");

%% Anonymous functions
oneil_calibration = @(d18O_ratio,seawater_d18O_ratio) ((((3.597122e-4)*log(d18O_ratio./seawater_d18O_ratio))+1.2194e-6).^(-0.5))-273.15;
kim_oneil_calibration = @(d18O_ratio,seawater_d18O_ratio) ((0.0554631*log(d18O_ratio./seawater_d18O_ratio) + 0.0017981).^(-1))-273.15;
no_ice_hansen_calibration = @(d18O) -4*d18O+12;
anderson_arthur_calibration = @(d18O,d18O_sw) 16-(4.14*(d18O-d18O_sw))+(0.13*(d18O-d18O_sw).^2);

%%
d18O = Geochemistry_Helpers.delta().create([height(data),1]);
d18O.assignToAll("standard","O_VPDB");
d18O.assignToEach("value",data.d18O);

seawater_d18O = Geochemistry_Helpers.delta("O_VSMOW",-1.2);

data.oneil_temperature = oneil_calibration(d18O.collate("ratio"),seawater_d18O.ratio);
data.kim_oneil_temperature = kim_oneil_calibration(d18O.collate("ratio"),seawater_d18O.ratio);
data.hansen_temperature = no_ice_hansen_calibration(data.d18O);
data.anderson_arthur_temperature = anderson_arthur_calibration(data.d18O,seawater_d18O.value);

%%
averaged.height = unique(data.height);
for depth_index = 1:numel(averaged.height)
    at_depth = data.height==averaged.height(depth_index);
    
    averaged.d13C(depth_index,1) = nanmean(data.d13C(at_depth));
    averaged.d18O(depth_index,1) = nanmean(data.d18O(at_depth));
    
    averaged.relative_age(depth_index,1) = nanmean(data.relative_age(at_depth));
    averaged.absolute_age(depth_index,1) = nanmean(data.absolute_age(at_depth));
    
    averaged.oneil_temperature(depth_index,1) = nanmean(data.oneil_temperature(at_depth));
    averaged.oneil_temperature_uncertainty(depth_index,1) = nanstd(data.oneil_temperature(at_depth));
    
    averaged.kim_oneil_temperature(depth_index,1) = nanmean(data.kim_oneil_temperature(at_depth));
    averaged.kim_oneil_temperature_uncertainty(depth_index,1) = nanstd(data.kim_oneil_temperature(at_depth));
    
    averaged.hansen_temperature(depth_index,1) = nanmean(data.hansen_temperature(at_depth));
    averaged.hansen_temperature_uncertainty(depth_index,1) = nanstd(data.hansen_temperature(at_depth));

    averaged.anderson_arthur_temperature(depth_index,1) = nanmean(data.anderson_arthur_temperature(at_depth));
    averaged.anderson_arthur_temperature_uncertainty(depth_index,1) = nanstd(data.anderson_arthur_temperature(at_depth));
end
averaged.oneil_temperature_uncertainty(averaged.oneil_temperature_uncertainty==0) = nanmean(averaged.oneil_temperature_uncertainty(averaged.oneil_temperature_uncertainty~=0));
averaged.kim_oneil_temperature_uncertainty(averaged.kim_oneil_temperature_uncertainty==0) = nanmean(averaged.kim_oneil_temperature_uncertainty(averaged.kim_oneil_temperature_uncertainty~=0));
averaged.hansen_temperature_uncertainty(averaged.hansen_temperature_uncertainty==0) = nanmean(averaged.hansen_temperature_uncertainty(averaged.hansen_temperature_uncertainty~=0));
averaged.anderson_arthur_temperature_uncertainty(averaged.anderson_arthur_temperature_uncertainty==0) = nanmean(averaged.anderson_arthur_temperature_uncertainty(averaged.anderson_arthur_temperature_uncertainty~=0));

averaged_table = struct2table(averaged);

%% Save
writematrix(["oneil_temperature","kim_oneil_temperature","hansen_temperature","anderson_arthur_temperature"],"./../Data/d18O_d13C.xlsx","Sheet","Matlab","Range","F1");
writematrix(["°C","°C","°C"],"./../Data/d18O_d13C.xlsx","Sheet","Matlab","Range","F2");
writematrix(data{:,:},"./../Data/d18O_d13C.xlsx","Sheet","Matlab","Range","A3");

writematrix(string(averaged_table.Properties.VariableNames),"./../Data/d18O_d13C.xlsx","Sheet","Averaged","Range","A1");
writematrix(["m","‰","‰","Myr","Ma","°C","°C","°C","°C","°C","°C","°C","°C"],"./../Data/d18O_d13C.xlsx","Sheet","Averaged","Range","A2");
writematrix(averaged_table{:,:},"./../Data/d18O_d13C.xlsx","Sheet","Averaged","Range","A3");
