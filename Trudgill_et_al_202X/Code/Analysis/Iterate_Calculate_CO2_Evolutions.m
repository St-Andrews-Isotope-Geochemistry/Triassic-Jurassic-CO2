% Run TJ CO2
data_directory = "./../../Data/";
filename = "TJ_CO2_Evolutions.csv";
for loop_index = 1:10
    clear all
    run("TJ_CO2.m");
    
    writematrix(pH,data_directory+filename,'WriteMode','append');
    writematrix(co2,data_directory+filename,'WriteMode','append');
    writematrix(saturation_state,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("dic"),data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("alkalinity"),data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("temperature"),data_directory+filename,'WriteMode','append');
    writematrix(d11B_out,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("calcium"),data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("magnesium"),data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.collate("boron").collate("epsilon"),data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.collate("boron").collate("d11B_sw").collate("value"),data_directory+filename,'WriteMode','append');
end