clear

%%
data_directory = "./../../../Data/";

%% Background d11B + Minimum Change
minimum_pH_change_inputs = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Input.json"));

%% Background pH + Minimum Change
minimum_pH_ensemble = readmatrix(data_directory+"/Minimum_pH_Change/TJ_Minimum_pH_Ensemble.csv");
pH = minimum_pH_ensemble(3:4,:);

background.pH = pH(1,:);
perturbation.pH = pH(2,:);

% If a value in either array is imaginary then we don't want it
real_pH = imag(perturbation.pH)==0 & imag(background.pH)==0;

assert(sum(real_pH)==numel(background.pH),"Some non real values.");

change.pH = perturbation.pH(real_pH) - background.pH(real_pH);
change.pH_distribution = Geochemistry_Helpers.Distribution.fromSamples(-1:0.01:1,change.pH).normalise();

background.pH_distribution = Geochemistry_Helpers.Distribution.fromSamples(7.5:0.01:9,background.pH).normalise();
background.pH_median = background.pH_distribution.median();
background.pH_uncertainty = background.pH_distribution.standard_deviation()*2;

change.minimum_pH = change.pH_distribution.quantile(0.95);
change.minimum_pH_uncertainty = change.pH_distribution.standard_deviation();

%% d11Bsw
d11B_sw = minimum_pH_ensemble(end,:);
d11B_sw_distribution = Geochemistry_Helpers.Distribution.fromSamples(20:0.1:40,d11B_sw).normalise();
d11B_sw_median = d11B_sw_distribution.median();
d11B_sw_uncertainty = d11B_sw_distribution.standard_deviation()*2;

%%
if exist(data_directory+"/Minimum_pH_Change/Manuscript_Values.json",'file')
    existing_manuscript_values = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Manuscript_Values.json"));
end
existing_manuscript_values.background_pH = round([background.pH_median,background.pH_uncertainty],3);
existing_manuscript_values.background_d11B = minimum_pH_change_inputs.initial_d11B;

existing_manuscript_values.minimum_pH_change = round([change.minimum_pH,change.minimum_pH_uncertainty],3);
existing_manuscript_values.minimum_d11B_change = minimum_pH_change_inputs.d11B_change;

existing_manuscript_values.d11B_sw = round([d11B_sw_median,d11B_sw_uncertainty],3);

writeInputJSON(data_directory+"/Minimum_pH_Change/Manuscript_Values.json",existing_manuscript_values);

% fileID = fopen(data_directory+"/Minimum_pH_Change/Manuscript_Values.json","w");
% fwrite(fileID,jsonencode(existing_manuscript_values));
% fclose(fileID);
