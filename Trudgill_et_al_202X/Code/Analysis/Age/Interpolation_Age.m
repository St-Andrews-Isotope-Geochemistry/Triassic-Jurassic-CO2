% Calculate ages to interpolate onto
clear

data_directory = "./../../../Data/";

boron_data = readtable(data_directory+"/Boron/TJ_d11B.xlsx","Sheet","With_Age");

interpolation_ages = unique(sort([boron_data.age',linspace(min(boron_data.age),max(boron_data.age),80)]));
structure.interpolation_ages = round(interpolation_ages,4);

%%
file = data_directory+"/Age/Interpolation_Age.json";
fileID = fopen(file,"w");
fwrite(fileID,jsonencode(structure));
fclose(fileID);
