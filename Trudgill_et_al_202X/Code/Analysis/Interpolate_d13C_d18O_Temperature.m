clear
data_directory = "./../../Data/";
d18O_d13C = readtable(data_directory+"TJ_d18O_d13C.xlsx","Sheet","Averaged");
d18O_d13C_delta = readtable(data_directory+"TJ_d18O_d13C.xlsx","Sheet","Delta_Temperature");
d11B = readtable(data_directory+"TJ_d11B.xlsx","Sheet","With_Age");
d11B_delta = readtable(data_directory+"TJ_d11B.xlsx","Sheet","With_Age");

%% Interpolate to same place as boron samples (linear)
d11B.d18O = interp1(d18O_d13C_delta.height,d18O_d13C_delta.d18O,d11B.height);
d11B.d18O_uncertainty = interp1(d18O_d13C_delta.height,d18O_d13C_delta.d18O_uncertainty,d11B.height);

d11B.d13C = interp1(d18O_d13C_delta.height,d18O_d13C_delta.d13C,d11B.height);
d11B.d13C_uncertainty = interp1(d18O_d13C_delta.height,d18O_d13C_delta.d13C_uncertainty,d11B.height);

d11B.oneil_temperature = interp1(d18O_d13C.height,d18O_d13C.oneil_temperature,d11B.height);
d11B.oneil_temperature_uncertainty = interp1(d18O_d13C.height,d18O_d13C.oneil_temperature_uncertainty,d11B.height);

d11B.kim_oneil_temperature = interp1(d18O_d13C.height,d18O_d13C.kim_oneil_temperature,d11B.height);
d11B.kim_oneil_temperature_uncertainty = interp1(d18O_d13C.height,d18O_d13C.kim_oneil_temperature_uncertainty,d11B.height);

d11B.hansen_temperature = interp1(d18O_d13C.height,d18O_d13C.hansen_temperature,d11B.height);
d11B.hansen_temperature_uncertainty = interp1(d18O_d13C.height,d18O_d13C.hansen_temperature_uncertainty,d11B.height);

d11B.anderson_arthur_temperature = interp1(d18O_d13C.height,d18O_d13C.anderson_arthur_temperature,d11B.height);
d11B.anderson_arthur_temperature_uncertainty = interp1(d18O_d13C.height,d18O_d13C.anderson_arthur_temperature_uncertainty,d11B.height);

%% Delta temperature
d11B_delta.d18O = d11B.d18O;
d11B_delta.d18O_uncertainty = d11B.d18O_uncertainty;

d11B_delta.d13C = d11B.d13C;
d11B_delta.d13C_uncertainty = d11B.d13C_uncertainty;

d11B_delta.delta_temperature = interp1(d18O_d13C_delta.height,d18O_d13C_delta.delta_temperature,d11B_delta.height);
d11B_delta.delta_temperature_uncertainty = interp1(d18O_d13C_delta.height,d18O_d13C_delta.delta_temperature_uncertainty,d11B_delta.height);

%% Interpolate to same place as boron samples (Gaussian)

%% Save
current_file = data_directory+"TJ_d11B.xlsx";
writematrix(string(d11B(:,2:end).Properties.VariableNames),current_file,"Sheet","Temperature_Calibrations","Range","A1");
writematrix(["m","Ma","‰","‰","‰","‰","‰","‰","°C","°C","°C","°C","°C","°C","°C","°C"],current_file,"Sheet","Temperature_Calibrations","Range","A2");
writematrix(d11B{:,2:end},current_file,"Sheet","Temperature_Calibrations","Range","A3");

%% Save
current_file = data_directory+"TJ_d11B.xlsx";
writematrix(string(d11B_delta(:,2:end).Properties.VariableNames),current_file,"Sheet","Delta_Temperature","Range","A1");
writematrix(["m","Ma","‰","‰","‰","‰","‰","‰","°C","°C"],current_file,"Sheet","Delta_Temperature","Range","A2");
writematrix(string(cell2mat(d11B_delta.sample)),current_file,"Sheet","Delta_Temperature","Range","A3");
writematrix(d11B_delta{:,2:end},current_file,"Sheet","Delta_Temperature","Range","A3");
