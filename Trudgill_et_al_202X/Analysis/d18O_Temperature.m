clear
d18O_d13C = readtable("./../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab");
d11B = readtable("./../Data/TJ_d11B_pH.xlsx","Sheet","Matlab");

%% Anonymous functions
oneil_calibration = @(d18O_ratio,seawater_d18O_ratio) ((((3.597122e-4)*log(d18O_ratio./seawater_d18O_ratio))+1.2194e-6).^(-0.5))-273.15;
kim_oneil_calibration = @(d18O_ratio,seawater_d18O_ratio) ((0.0554631*log(d18O_ratio./seawater_d18O_ratio) + 0.0017981).^(-1))-273.15;
no_ice_hansen_calibration = @(d18O) -4*d18O+12;
anderson_arthur_calibration = @(d18O,d18O_sw) 16-(4.14*(d18O-d18O_sw))+(0.13*(d18O-d18O_sw).^2);

%%
d18O = Geochemistry_Helpers.delta().create([height(d18O_d13C),1]);
d18O.assignToAll("standard","O_VPDB");
d18O.assignToEach("value",d18O_d13C.d18O);

seawater_d18O = Geochemistry_Helpers.delta("O_VSMOW",-1.2);

d18O_d13C.oneil_temperature = oneil_calibration(d18O.collate("ratio"),seawater_d18O.ratio);
d18O_d13C.kim_oneil_temperature = kim_oneil_calibration(d18O.collate("ratio"),seawater_d18O.ratio);
d18O_d13C.hansen_temperature = no_ice_hansen_calibration(d18O_d13C.d18O);
d18O_d13C.anderson_arthur_temperature = anderson_arthur_calibration(d18O_d13C.d18O,seawater_d18O.value);

%%
averaged.height = unique(d18O_d13C.height);
for depth_index = 1:numel(averaged.height)
    at_depth = d18O_d13C.height==averaged.height(depth_index);
    
    averaged.d13C(depth_index,1) = nanmean(d18O_d13C.d13C(at_depth));
    averaged.d18O(depth_index,1) = nanmean(d18O_d13C.d18O(at_depth));
    
    averaged.relative_age(depth_index,1) = nanmean(d18O_d13C.relative_age(at_depth));
    averaged.absolute_age(depth_index,1) = nanmean(d18O_d13C.absolute_age(at_depth));
    
    averaged.oneil_temperature(depth_index,1) = nanmean(d18O_d13C.oneil_temperature(at_depth));
    averaged.oneil_temperature_uncertainty(depth_index,1) = nanstd(d18O_d13C.oneil_temperature(at_depth));
    
    averaged.kim_oneil_temperature(depth_index,1) = nanmean(d18O_d13C.kim_oneil_temperature(at_depth));
    averaged.kim_oneil_temperature_uncertainty(depth_index,1) = nanstd(d18O_d13C.kim_oneil_temperature(at_depth));
    
    averaged.hansen_temperature(depth_index,1) = nanmean(d18O_d13C.hansen_temperature(at_depth));
    averaged.hansen_temperature_uncertainty(depth_index,1) = nanstd(d18O_d13C.hansen_temperature(at_depth));

    averaged.anderson_arthur_temperature(depth_index,1) = nanmean(d18O_d13C.anderson_arthur_temperature(at_depth));
    averaged.anderson_arthur_temperature_uncertainty(depth_index,1) = nanstd(d18O_d13C.anderson_arthur_temperature(at_depth));
end
averaged.oneil_temperature_uncertainty(averaged.oneil_temperature_uncertainty==0) = nanmean(averaged.oneil_temperature_uncertainty(averaged.oneil_temperature_uncertainty~=0));
averaged.kim_oneil_temperature_uncertainty(averaged.kim_oneil_temperature_uncertainty==0) = nanmean(averaged.kim_oneil_temperature_uncertainty(averaged.kim_oneil_temperature_uncertainty~=0));
averaged.hansen_temperature_uncertainty(averaged.hansen_temperature_uncertainty==0) = nanmean(averaged.hansen_temperature_uncertainty(averaged.hansen_temperature_uncertainty~=0));
averaged.anderson_arthur_temperature_uncertainty(averaged.anderson_arthur_temperature_uncertainty==0) = nanmean(averaged.anderson_arthur_temperature_uncertainty(averaged.anderson_arthur_temperature_uncertainty~=0));

averaged = struct2table(averaged);

%% Interpolate to same place as boron samples
d11B.oneil_temperature = interp1(averaged.height,averaged.oneil_temperature,d11B.height);
d11B.oneil_temperature_uncertainty = interp1(averaged.height,averaged.oneil_temperature_uncertainty,d11B.height);

d11B.kim_oneil_temperature = interp1(averaged.height,averaged.kim_oneil_temperature,d11B.height);
d11B.kim_oneil_temperature_uncertainty = interp1(averaged.height,averaged.kim_oneil_temperature_uncertainty,d11B.height);

d11B.hansen_temperature = interp1(averaged.height,averaged.hansen_temperature,d11B.height);
d11B.hansen_temperature_uncertainty = interp1(averaged.height,averaged.hansen_temperature_uncertainty,d11B.height);

d11B.anderson_arthur_temperature = interp1(averaged.height,averaged.anderson_arthur_temperature,d11B.height);
d11B.anderson_arthur_temperature_uncertainty = interp1(averaged.height,averaged.anderson_arthur_temperature_uncertainty,d11B.height);

%% Maximum initial temperature
maximum_initial_temperature = nanmax(d11B.hansen_temperature(1:9));
maximum_initial_temperature_row = d11B(d11B.hansen_temperature(1:9)==maximum_initial_temperature,:);
maximum_initial_temperatures = [maximum_initial_temperature_row.oneil_temperature,maximum_initial_temperature_row.kim_oneil_temperature,maximum_initial_temperature_row.hansen_temperature,maximum_initial_temperature_row.anderson_arthur_temperature];
maximum_initial_temperature_centre = mean(maximum_initial_temperatures);
maximum_initial_temperature_uncertainty = max(maximum_initial_temperatures)-min(maximum_initial_temperatures);

minimum_temperature_change = maximum_initial_temperature-d11B.hansen_temperature(d11B.d11B(10:end)==nanmin(d11B.d11B(10:end)));
minimum_temperature_change_uncertainty = d11B.hansen_temperature_uncertainty(d11B.d11B(10:end)==nanmin(d11B.d11B(10:end)));

%% Save
writematrix(["oneil_temperature","kim_oneil_temperature","hansen_temperature","anderson_arthur_temperature"],"./../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab","Range","F1");
writematrix(["°C","°C","°C","°C"],"./../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab","Range","F2");
writematrix(d18O_d13C{:,end-3:end},"./../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab","Range","F3");

writematrix(string(averaged.Properties.VariableNames),"./../Data/TJ_d18O_d13C.xlsx","Sheet","Averaged","Range","A1");
writematrix(["m","‰","‰","Myr","Ma","°C","°C","°C","°C","°C","°C","°C","°C"],"./../Data/TJ_d18O_d13C.xlsx","Sheet","Averaged","Range","A2");
writematrix(averaged{:,:},"./../Data/TJ_d18O_d13C.xlsx","Sheet","Averaged","Range","A3");

writematrix(["oneil_temperature","oneil_temperature_uncertainty","kim_oneil_temperature","kim_oneil_temperature_uncertainty","hansen_temperature","hansen_temperature_uncertainty","anderson_arthur_temperature","anderson_arthur_temperature_uncertainty"],"./../Data/TJ_d11B_pH.xlsx","Sheet","Matlab","Range","G1");
writematrix(["°C","°C","°C","°C","°C","°C","°C","°C"],"./../Data/TJ_d11B_pH.xlsx","Sheet","Matlab","Range","G2");
writematrix(d11B{:,end-7:end},"./../Data/TJ_d11B_pH.xlsx","Sheet","Matlab","Range","G3");

maximum_initial_temperature_file = fopen("./../Data/Maximum_Initial_Temperature.txt","w");
fprintf(maximum_initial_temperature_file,"Maximum initial temperature = "+num2str(maximum_initial_temperature_centre));
fprintf(maximum_initial_temperature_file," +- "+num2str(maximum_initial_temperature_uncertainty));
fprintf(maximum_initial_temperature_file,char(176)+"C"+newline);
fprintf(maximum_initial_temperature_file,"Minimum temperature change = +"+num2str(minimum_temperature_change));
fprintf(maximum_initial_temperature_file," +- "+num2str(minimum_temperature_change_uncertainty));
fprintf(maximum_initial_temperature_file,char(176)+"C");
fclose(maximum_initial_temperature_file);
