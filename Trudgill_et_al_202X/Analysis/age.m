clear
data = readtable("./d18O_d13C","Sheet","Matlab");
age_data = readtable("./d18O_d13C","Sheet","Age_Map");

ages = interp1(age_data.Height_Korte,age_data.Relative_Age_Korte,data.Depth)*1000;
absolute_ages = (201486.2222222222-ages)/1e3;

