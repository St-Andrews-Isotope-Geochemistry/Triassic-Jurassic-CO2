clear

%% Load in the data
d18O_d13C = readtable("./../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab");
d11B = readtable("./../Data/TJ_d11B_pH.xlsx","Sheet","Matlab");

age_data = readtable("./../Data/TJ_Age_Calibration.xlsx","Sheet","Matlab");

%% Interpolate ages based on tiepoints
d18O_d13C.relative_age = interp1(age_data.height_CVR,age_data.age_CIE,d18O_d13C.height)/1e3;
d18O_d13C.absolute_age = 201.564-d18O_d13C.relative_age;

d11B.relative_age = interp1(age_data.height_CVR,age_data.age_CIE,d11B.height)/1e3;
d11B.absolute_age = 201.564-d11B.relative_age;

%% Save
writematrix(["relative_age","absolute_age"],"./../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab","Range","D1");
writematrix(["Myr","Ma"],"./../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab","Range","D2");
writematrix(d18O_d13C{:,end-1:end},"./../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab","Range","D3");

writematrix(["relative_age","absolute_age"],"./../Data/TJ_d11B_pH.xlsx","Sheet","Matlab","Range","E1");
writematrix(["Myr","Ma"],"./../Data/TJ_d11B_pH.xlsx","Sheet","Matlab","Range","E2");
writematrix(d11B{:,end-1:end},"./../Data/TJ_d11B_pH.xlsx","Sheet","Matlab","Range","E3");

