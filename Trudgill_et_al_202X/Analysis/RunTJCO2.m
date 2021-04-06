% Run TJ CO2

for loop_index = 1:8
    clear all
    run("TJ_CO2.m");
    
    writematrix(pH,"Output.csv",'WriteMode','append');
    writematrix(co2,"Output.csv",'WriteMode','append');
    writematrix(saturation_state,"Output.csv",'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("dic"),"Output.csv",'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("alkalinity"),"Output.csv",'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("temperature"),"Output.csv",'WriteMode','append');
    writematrix(d11B_out,"Output.csv",'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("calcium"),"Output.csv",'WriteMode','append');
    writematrix(co2_evolutions.collate("carbonate_chemistry").collate("magnesium"),"Output.csv",'WriteMode','append');
    writematrix(co2_evolutions.collate("boron").collate("epsilon"),"Output.csv",'WriteMode','append');
    writematrix(co2_evolutions.collate("boron").collate("d11B_sw").collate("value"),"Output.csv",'WriteMode','append');
end