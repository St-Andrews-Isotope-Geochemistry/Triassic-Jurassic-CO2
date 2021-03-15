clear

%% Load in the data
data = readtable("./../Data/d18O_d13C","Sheet","Matlab");
age_data = readtable("./../Data/d18O_d13C","Sheet","Age_Calibration");

%% Interpolate ages based on tiepoints
data.relative_age = interp1(age_data.height_CVR,age_data.age_CIE,data.height)/1e3;
data.absolute_age = 201.564-data.relative_age;

%% Save
writematrix(["relative_age","absolute_age"],"./../Data/d18O_d13C.xlsx","Sheet","Matlab","Range","D1");
writematrix(["Myr","Ma"],"./../Data/d18O_d13C.xlsx","Sheet","Matlab","Range","D2");
writematrix(data{:,:},"./../Data/d18O_d13C.xlsx","Sheet","Matlab","Range","A3");

