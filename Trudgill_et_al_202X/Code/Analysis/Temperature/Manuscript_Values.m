clear

%%
data_directory = "./../../../Data/";

%% Delta temperature
minimum_pH_change_inputs = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Input.json"));

%%
if exist(data_directory+"/Temperature/Manuscript_Values.json",'file')
    existing_manuscript_values = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Manuscript_Values.json"));
end
existing_manuscript_values.temperature_change = minimum_pH_change_inputs.temperature_change;

writeInputJSON(data_directory+"/Temperature/Manuscript_Values.json",existing_manuscript_values);

