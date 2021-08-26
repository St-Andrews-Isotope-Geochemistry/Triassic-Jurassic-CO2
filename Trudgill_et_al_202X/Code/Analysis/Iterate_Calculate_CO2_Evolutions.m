% Run TJ CO2
clear
data_directory = "./../../Data/";
filename = "TJ_CO2_Evolutions.csv";

delete(data_directory+filename);
    
for loop_index = 1:10
    clearvars -except loop_index data_directory filename    

    run("Calculate_CO2_Evolutions.m");    

    writematrix(pH_output,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions_value,data_directory+filename,'WriteMode','append');
    writematrix(saturation_state_output,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.carbonate_chemistry.dic,data_directory+filename,'WriteMode','append');
    writematrix(alkalinity_output,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.carbonate_chemistry.temperature,data_directory+filename,'WriteMode','append');
    writematrix(d11B_evolutions_viable,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.carbonate_chemistry.calcium,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.carbonate_chemistry.magnesium,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.boron.epsilon,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.boron.d11B_sw.collate("value"),data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.species_calibration.coefficients(:,:,1),data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.species_calibration.coefficients(:,:,2),data_directory+filename,'WriteMode','append');
end