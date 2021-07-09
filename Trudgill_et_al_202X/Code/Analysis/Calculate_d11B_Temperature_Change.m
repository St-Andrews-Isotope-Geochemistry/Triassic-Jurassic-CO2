%% Calculate d11B and temperature change
% To calculate the minimum ΔpH we need to know dependent properties:
% - Preperturbation d11B
% - Change in measured d11B
% - Preperturbation temperature
% - Change in temperature
%
% Preperturbation
% For both d11B and temperature we take the average of the first 9
% datapoints (for temperature this has been interpolated from the original)
%
% Perturbation
% For both d11B and temperature we take the average of the 2 datapoints
% with the lowest d11B and their associated temperature

clear
data_directory = "./../../Data/";
d11B = readtable(data_directory+"TJ_d11B.xlsx","Sheet","Temperature_Calibrations");
d11B_delta = readtable(data_directory+"TJ_d11B.xlsx","Sheet","Delta_Temperature");

d11B_combined = [d11B,d11B_delta(:,end-1:end)];

%% Separate data
d11B_preperturbation = d11B_combined(1:9,:);
d11B_sorted = sortrows(d11B_combined,"d11B","ascend");
d11B_minima = d11B_sorted(1:2,:);

%% Calculate metrics
% d11B
preperturbation.d11B = mean(d11B_preperturbation.d11B);
preperturbation.d11B_uncertainty = sqrt(sum(d11B_preperturbation.d11B_uncertainty.^2))/height(d11B_preperturbation);

perturbation.d11B = mean(d11B_minima.d11B);
perturbation.d11B_uncertainty = sqrt(sum(d11B_minima.d11B_uncertainty.^2))/height(d11B_minima);

% Temperature
preperturbation.temperature = mean([d11B_preperturbation.oneil_temperature,d11B_preperturbation.kim_oneil_temperature,d11B_preperturbation.hansen_temperature,d11B_preperturbation.anderson_arthur_temperature]);
preperturbation.temperature_uncertainty = sqrt(sum([d11B_preperturbation.oneil_temperature_uncertainty,d11B_preperturbation.kim_oneil_temperature_uncertainty,d11B_preperturbation.hansen_temperature_uncertainty,d11B_preperturbation.anderson_arthur_temperature_uncertainty].^2))/height(d11B_preperturbation);

% We have four calibrations so having calculated the metrics for each, we
% combine into a single probability distribution (assuming each calibration
% is equally likely to be correct)
for index = 1:numel(preperturbation.temperature)
    temperature_distributions(index) = Geochemistry_Helpers.Distribution(0:0.1:30,"Gaussian",[preperturbation.temperature(index),preperturbation.temperature_uncertainty(index)]).normalise();
end
initial_temperature = Geochemistry_Helpers.Distribution(0:0.1:30,"Manual",sum(temperature_distributions.probabilities)).normalise();

% The initial temperature distribution is bimodal. We model that with a
% combination of two Gaussians.
fit_function = @(a1,b1,c1,a2,b2,c2,x) a1*exp(-((x-b1)/c1).^2) + a2*exp(-((x-b2)/c2).^2);
initial_temperature_fit = fit(initial_temperature.bin_midpoints,initial_temperature.probabilities,fit_function,'StartPoint',[0.03,10,1,0.03,13,1],'Robust','LAR');

% For the change + uncertainty we use the offset from the first datapoint
temperature_change = mean(d11B_minima.delta_temperature);
temperature_change_uncertainty = sqrt(sum(d11B_minima.delta_temperature_uncertainty.^2))/height(d11B_minima);

%% Separate into named parameters
initial_d11B = preperturbation.d11B;
initial_d11B_uncertainty = preperturbation.d11B_uncertainty;

maximum_initial_d11B = initial_d11B+(2*initial_d11B_uncertainty);
maximum_initial_d11B_uncertainty = initial_d11B_uncertainty/2;

d11B_change = perturbation.d11B-preperturbation.d11B;
d11B_change_uncertainty = sqrt(sum([preperturbation.d11B_uncertainty,perturbation.d11B_uncertainty].^2));

minimum_d11B_change = d11B_change-(2*d11B_change_uncertainty);
minimum_d11B_change_uncertainty = d11B_change_uncertainty/2;

initial_temperature_parameters = reshape(coeffvalues(initial_temperature_fit),3,2)';
maximum_initial_temperature = initial_temperature.quantile(0.95);
maximum_initial_temperature_uncertainty = initial_temperature.standard_deviation/2;

minimum_temperature_change = temperature_change-(2*temperature_change_uncertainty);
minimum_temperature_change_uncertainty = temperature_change_uncertainty./2;

%% Save
current_file_contents = jsondecode(fileread(data_directory+"Minimum_pH_Change_Input.json"));
current_file_contents.initial_d11B = [maximum_initial_d11B,maximum_initial_d11B_uncertainty];
current_file_contents.d11B_change = [minimum_d11B_change,minimum_d11B_change_uncertainty];

current_file_contents.initial_temperature = [maximum_initial_temperature,maximum_initial_temperature_uncertainty];
current_file_contents.temperature_change = [minimum_temperature_change,minimum_temperature_change_uncertainty];

current_file_fieldnames = string(fieldnames(current_file_contents));

fileID = fopen(data_directory+"Minimum_pH_Change_Input.json","w");
fwrite(fileID,"{"+newline);
for name = current_file_fieldnames'
    fwrite(fileID,string(char(9))+'"'+name+'":');
    if numel(current_file_contents.(name))==1
        fwrite(fileID,num2str(current_file_contents.(name)));
    elseif numel(current_file_contents.(name))==2
        values = current_file_contents.(name);
        fwrite(fileID,"["+num2str(values(1),3)+","+num2str(values(2),3)+"]");        
    end
    
    if ~strcmp(name,current_file_fieldnames(end))
        fwrite(fileID,","+newline);
    else
        fwrite(fileID,newline+"}");
    end
end
fclose(fileID);
