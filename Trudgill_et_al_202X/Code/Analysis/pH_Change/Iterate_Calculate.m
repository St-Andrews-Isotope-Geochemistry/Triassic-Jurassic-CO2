% Run TJ CO2
clear
data_directory = "./../../../Data/";
interpolation_ages = jsondecode(fileread(data_directory+"/Age/Interpolation_Age.json")).interpolation_ages;

createResultsFile(data_directory+"/pH_Change/prior.json","prior",interpolation_ages');
createResultsFile(data_directory+"/pH_Change/intermediate.json","intermediate",interpolation_ages');
createResultsFile(data_directory+"/pH_Change/posterior.json","posterior",interpolation_ages');

for loop_index = 1
    clearvars -except loop_index data_directory filename    
    
    run("Calculate.m");
    
    tab = string(char(9));
    for prior_index = 1:round_1.number_of_samples
        insertSample(data_directory+"/pH_Change/prior.json",round_1,prior_index);
    end
    
    for intermediate_index = 1:round_2.number_of_samples
        insertSample(data_directory+"/pH_Change/intermediate.json",round_2,intermediate_index);
    end
    
    for posterior_index = 1:round_3.number_of_samples
        insertSample(data_directory+"/pH_Change/posterior.json",round_3,posterior_index);
        insertpHTimeSeries(data_directory+"/ph_Change/posterior.json",round_3,posterior_index);        
    end
    
    fclose("all");
end

%%
function insertSample(file,round,index)
    tab = string(char(9));
    fileID = fopen(file,"r+");
    fseek(fileID,-6,'eof');
    if strcmp(fgetl(fileID),"[")
        fwrite(fileID,tab+tab+"{"+tab+tab);
    else
        fseek(fileID,-5,'eof');
        fwrite(fileID,","+newline+tab+tab+"{"+tab+tab);
    end
    insertValue(fileID,"co2_initial",round.co2.initial_sampler.samples(index)*1e6);
    insertValue(fileID,"saturation_state_initial",round.saturation_state.initial_sampler.samples(index));
    insertValue(fileID,"pH_initial",round.pH.initial_samples(index));
    insertValue(fileID,"calcium",round.calcium.sampler.samples(index));
    insertValue(fileID,"magnesium",round.magnesium.sampler.samples(index));
    insertValue(fileID,"epsilon",round.epsilon.sampler.samples(index));
    insertArray(fileID,"species_calibration",[round.species_calibration.gradient_sampler.samples(index),round.species_calibration.intercept_sampler.samples(index)]);
    insertValue(fileID,"d11B_sw",round.d11B_sw.sampler.samples(index));
    insertArray(fileID,"d11B_measured",round.d11B_measured.samples(:,index)');
    insertArray(fileID,"d11B_4",round.d11B_4.samples(:,index)');
    insertArray(fileID,"temperature",round.temperature.samples(:,index)');
    insertArray(fileID,"alkalinity",round.alkalinity.samples(:,index)'*1e6);
    fseek(fileID,-1,'eof');
    
    fwrite(fileID,newline+tab+tab+"}"+newline+tab+"]"+newline+"}");
    fclose(fileID);
end
function createResultsFile(file,type,interpolation_ages)
    tab = string(char(9));
    fileID = fopen(file,"w");
    fwrite(fileID,"{"+newline);
    fwrite(fileID,tab+'"'+"age"+'"'+" : ["+strip(strrep(num2str(interpolation_ages,"%.5G,")," ",""),',')+"],");
    fwrite(fileID,newline+tab+'"'+type+'"'+": ["+newline+tab+"]"+newline+"}");
    fclose(fileID);
end
function insertValue(fileID,name,value)
    tab = string(char(9));
    fwrite(fileID,newline+tab+tab+tab);
    fwrite(fileID,'"'+name+'"'+" : "+num2str(value,"%.5G,"));
end
function insertArray(fileID,name,value)
    tab = string(char(9));
    fwrite(fileID,newline+tab+tab+tab);
    fwrite(fileID,'"'+name+'"'+" : ["+strip(strrep(num2str(value,"%.5G,")," ",""),',')+"],");
end
function insertpHTimeSeries(file,round,index)
    tab = string(char(9));
    fileID = fopen(file,"r+");
    fseek(fileID,-5,'eof');
    fwrite(fileID,","+newline+tab+tab+"{"+tab+tab);
    insertArray(fileID,"pH",round.pH.samples(:,index)');
    fclose(fileID);
end