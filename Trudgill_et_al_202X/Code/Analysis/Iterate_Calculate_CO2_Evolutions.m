% Run TJ CO2
for loop_index = 1:1
    clear all
    run("Calculate_CO2_Evolutions.m");
    
    data_directory = "./../../Data/";
    filename = "TJ_CO2_Evolutions.csv";
    
    delete(data_directory+filename);

    writematrix(pH,data_directory+filename,'WriteMode','append');
    writematrix(co2,data_directory+filename,'WriteMode','append');
    writematrix(saturation_state,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.carbonate_chemistry.dic,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.carbonate_chemistry.alkalinity,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.carbonate_chemistry.temperature,data_directory+filename,'WriteMode','append');
    writematrix(d11B_out,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.carbonate_chemistry.calcium,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.carbonate_chemistry.magnesium,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.boron.epsilon,data_directory+filename,'WriteMode','append');
    writematrix(co2_evolutions.boron.d11B_sw.collate("value"),data_directory+filename,'WriteMode','append');
end