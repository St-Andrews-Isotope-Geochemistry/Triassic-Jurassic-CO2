clear

%% Load in the data
data_directory = "./../../../Data/";

d18O_d13C = readtable(data_directory+"/Temperature/TJ_d18O_d13C.xlsx","Sheet","Reformatted");
d11B = readtable(data_directory+"/Boron/TJ_d11B.xlsx","Sheet","Reformatted");

age_correlation_data = readtable(data_directory+"/Age/TJ_Age_Calibration.xlsx","Sheet","Correlation");
age_cyclostratigraphy_data = readtable(data_directory+"/Age/TJ_Age_Calibration.xlsx","Sheet","Cyclostratigraphy");

%% Interpolate ages based on tiepoints
age_cyclostratigraphy_data.height_CVR = interp1(age_correlation_data.height_hesselbo,age_correlation_data.height_CVR,age_cyclostratigraphy_data.height);
age_cyclostratigraphy_data_clean = rmmissing(age_cyclostratigraphy_data);

d18O_d13C_relative_age = (interp1(age_cyclostratigraphy_data_clean.height_CVR,age_cyclostratigraphy_data_clean.age,d18O_d13C.height)+77.77777777)/1e3;
d18O_d13C.age = 201.564-d18O_d13C_relative_age;

d11B_relative_age = (interp1(age_cyclostratigraphy_data_clean.height_CVR,age_cyclostratigraphy_data_clean.age,d11B.height)+77.77777777)/1e3;
d11B.age = 201.564-d11B_relative_age;

%% Save
current_file = data_directory+"/Temperature/TJ_d18O_d13C.xlsx";
writematrix(["height","age","d13C","d18O"],current_file,"Sheet","With_Age","Range","A1");
writematrix(["m","Ma","‰","‰"],current_file,"Sheet","With_Age","Range","A2");
writematrix([d18O_d13C.height,d18O_d13C.age,d18O_d13C.d13C,d18O_d13C.d18O],current_file,"Sheet","With_Age","Range","A3");

current_file = data_directory+"/Boron/TJ_d11B.xlsx";
writematrix(["sample","height","age","d11B","d11B_uncertainty"],current_file,"Sheet","With_Age","Range","A1");
writematrix([" ","m","Ma","‰","‰"],current_file,"Sheet","With_Age","Range","A2");
writematrix(string(cell2mat(d11B.sample)),current_file,"Sheet","With_Age","Range","A3");
writematrix([d11B.height,d11B.age,d11B.d11B,d11B.d11B_uncertainty],current_file,"Sheet","With_Age","Range","B3");

