clear
data_directory = "./../../Data/";
d18O_d13C = readtable(data_directory+"TJ_d18O_d13C.xlsx","Sheet","Delta_Temperature");
d11B = readtable(data_directory+"TJ_d11B.xlsx","Sheet","With_Age");

%% Interpolate to same place as boron samples (linear)
d11B.d18O = interp1(d18O_d13C.height,d18O_d13C.d18O,d11B.height);
d11B.d18O_uncertainty = interp1(d18O_d13C.height,d18O_d13C.d18O_uncertainty,d11B.height);

d11B.d13C = interp1(d18O_d13C.height,d18O_d13C.d13C,d11B.height);
d11B.d13C_uncertainty = interp1(d18O_d13C.height,d18O_d13C.d13C_uncertainty,d11B.height);

d11B.delta_temperature = interp1(d18O_d13C.height,d18O_d13C.delta_temperature,d11B.height);
d11B.delta_temperature_uncertainty = interp1(d18O_d13C.height,d18O_d13C.delta_temperature_uncertainty,d11B.height);

%% Interpolate to same place as boron samples (Gaussian)


%% Save
current_file = data_directory+"TJ_d11B.xlsx";
writematrix(string(d11B.Properties.VariableNames),current_file,"Sheet","With_Temperature","Range","A1");
writematrix([" ","m","Ma","‰","‰","‰","‰","‰","‰","°C","°C"],current_file,"Sheet","With_Temperature","Range","A2");
writematrix(string(cell2mat(d11B.sample)),current_file,"Sheet","With_Temperature","Range","A3");
writematrix(d11B{:,2:end},current_file,"Sheet","With_Temperature","Range","B3");
