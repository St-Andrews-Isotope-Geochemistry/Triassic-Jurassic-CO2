clear

%% Load data
data_directory = "./../../../Data/";
boron_data = readtable(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx","Sheet","Temperature_Calibrations");

posterior = jsondecode(fileread(data_directory+"/pH_change/posterior.json"));

interpolation_ages = jsondecode(fileread(data_directory+"/Age/Interpolation_Age.json")).interpolation_ages;
age_counter = 0;
for age = round(linspace(min(boron_data.age),max(boron_data.age),80),4)
    matching_age = find(interpolation_ages-age==0);
    if ~isnan(matching_age)
        age_counter = age_counter+1;
        age_location(age_counter) = matching_age;
    end
end

raw_evolutions_metrics = readmatrix(data_directory+"/pH_Change/Metrics_standard_smooth.csv");

number_of_quantiles = size(raw_evolutions_metrics,1)/4;
d11B_measured_quantiles = raw_evolutions_metrics(1:number_of_quantiles,age_location)';
pH_quantiles = raw_evolutions_metrics(number_of_quantiles+1:2*number_of_quantiles,age_location)';
saturation_state_quantiles = raw_evolutions_metrics(2*number_of_quantiles+1:3*number_of_quantiles,age_location)';
co2_quantiles = raw_evolutions_metrics(3*number_of_quantiles+1:end,age_location)';

% raw_minimum_metrics = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Metrics.json"));
raw_minimum_metrics = readmatrix(data_directory+"/Minimum_pH_Change/Quantiles.csv");
raw_evolutions_metrics_delta = readmatrix(data_directory+"/pH_Change/Metrics_Delta_standard_smooth.csv");

minimum_pH_change_d11B_measured = raw_minimum_metrics(1:2,:);
minimum_pH_change_pH = raw_minimum_metrics(3:4,:);

delta_pH_quantiles = raw_evolutions_metrics_delta(1:number_of_quantiles,age_location)';
delta_saturation_state_quantiles = raw_evolutions_metrics_delta(number_of_quantiles+1:2*number_of_quantiles,age_location)';
delta_co2_quantiles = raw_evolutions_metrics_delta(2*number_of_quantiles+1:end,age_location)';

%% Make a plot like figure 2 in the paper
all_pH.colour = [0.6,0,0];
high_initial_pH.colour = [0,0.5,0.8];
medium_initial_co2.colour = [0.8,0.8,0.2];
low_initial_co2.colour = [0,0,0.6];
marker_size = 4;
line_width = 1.5;

initial_age = [boron_data.age(1),boron_data.age(9)];
perturbation_age = [201.32,201.29];

figure(1);
clf
% subplot_handles(1) = subplot(4,1,1);
% hold on
% for quantile_index = 1:floor(number_of_quantiles/2)
%     upper = d11B_measured_quantiles(:,100-quantile_index);
%     lower = d11B_measured_quantiles(:,quantile_index);
%     patch([interpolation_ages;flipud(interpolation_ages)],[lower;flipud(upper)],all_pH.colour,'FaceAlpha',0.026,'EdgeColor','None');
% end
% plot(interpolation_ages,d11B_measured_quantiles(:,ceil(number_of_quantiles/2)),'-','Color',all_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
% 
% for boron_index = 1:height(boron_data)
%     plot([boron_data.age(boron_index),boron_data.age(boron_index)],boron_data.d11B(boron_index)+2.*[-boron_data.d11B_uncertainty(boron_index),boron_data.d11B_uncertainty(boron_index)],'k-');
% end
% d11B_plot_handles(1) = plot(boron_data.age,boron_data.d11B,'k.');
% 
% for quantile_index = 1:floor(number_of_quantiles/2)
%     upper = minimum_pH_change_d11B_measured(1,100-quantile_index);
%     lower = minimum_pH_change_d11B_measured(1,quantile_index);
%     patch([initial_age,fliplr(initial_age)],[lower,lower,upper,upper],[0.5,0.5,0.5],'FaceAlpha',0.1,'EdgeColor','None');
%     
%     upper = minimum_pH_change_d11B_measured(2,100-quantile_index);
%     lower = minimum_pH_change_d11B_measured(2,quantile_index);
%     patch([perturbation_age,fliplr(perturbation_age)],[lower,lower,upper,upper],[0.5,0.5,0.5],'FaceAlpha',0.1,'EdgeColor','None');
% end
% plot(initial_age,[minimum_pH_change_d11B_measured(1,ceil(number_of_quantiles/2)),minimum_pH_change_d11B_measured(1,ceil(number_of_quantiles/2))],'Color','k','LineWidth',0.3);
% d11B_plot_handles(2) = plot(perturbation_age,[minimum_pH_change_d11B_measured(2,ceil(number_of_quantiles/2)),minimum_pH_change_d11B_measured(2,ceil(number_of_quantiles/2))],'Color','k','LineWidth',0.3);
% 
% 
% ylabel("\delta^{11}B_{measured}");
% 
% legend(d11B_plot_handles,["Original data","Minimum \Delta\delta^{11}B"],'Location','SouthWest');

subplot_handles(1) = subplot(3,1,1);
hold on
for quantile_index = 1:floor(number_of_quantiles/2)
    upper = pH_quantiles(:,100-quantile_index);
    lower = pH_quantiles(:,quantile_index);
    patch([interpolation_ages(age_location);flipud(interpolation_ages(age_location))],[lower;flipud(upper)],all_pH.colour,'FaceAlpha',0.026,'EdgeColor','None');
end
pH_handles(1) = plot(interpolation_ages(age_location),pH_quantiles(:,ceil(number_of_quantiles/2)),'-','Color',all_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);


for quantile_index = 1:floor(number_of_quantiles/2)
    upper = minimum_pH_change_pH(1,100-quantile_index);
    lower = minimum_pH_change_pH(1,quantile_index);
    patch([initial_age,fliplr(initial_age)],[lower,lower,upper,upper],[0.5,0.5,0.5],'FaceAlpha',0.1,'EdgeColor','None');
    
    upper = minimum_pH_change_pH(2,100-quantile_index);
    lower = minimum_pH_change_pH(2,quantile_index);
    patch([perturbation_age,fliplr(perturbation_age)],[lower,lower,upper,upper],[0.5,0.5,0.5],'FaceAlpha',0.1,'EdgeColor','None');
end
plot(initial_age,[minimum_pH_change_pH(1,ceil(number_of_quantiles/2)),minimum_pH_change_pH(1,ceil(number_of_quantiles/2))],'Color','k','LineWidth',1);
pH_handles(2) = plot(perturbation_age,[minimum_pH_change_pH(2,ceil(number_of_quantiles/2)),minimum_pH_change_pH(2,ceil(number_of_quantiles/2))],'Color','k','LineWidth',1);

% patch([initial_age,fliplr(initial_age)],[raw_minimum_metrics.initial.pH_025,raw_minimum_metrics.initial.pH_025,raw_minimum_metrics.initial.pH_975,raw_minimum_metrics.initial.pH_975],[0.7,0.7,0.7],'FaceAlpha',0.8,'EdgeColor','None');

% patch([perturbation_age,fliplr(perturbation_age)],[raw_minimum_metrics.after.pH_025,raw_minimum_metrics.after.pH_025,raw_minimum_metrics.after.pH_975,raw_minimum_metrics.after.pH_975],[0.7,0.7,0.7],'FaceAlpha',0.8,'EdgeColor','None');

ylabel("pH");

% legend_handle = legend(pH_handles,["All evolutions","Minimum \DeltapH"],'Location','SouthWest');

subplot_handles(2) = subplot(3,1,2);
hold on
for quantile_index = 1:floor(number_of_quantiles/2)
    upper = delta_co2_quantiles(:,100-quantile_index);
    lower = delta_co2_quantiles(:,quantile_index);
    patch([interpolation_ages(age_location);flipud(interpolation_ages(age_location))],[lower;flipud(upper)],all_pH.colour,'FaceAlpha',0.026,'EdgeColor','None');
end
plot(interpolation_ages(age_location),(delta_co2_quantiles(:,ceil(number_of_quantiles/2))),'-','Color',all_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);

ylabel("CO_2 (ppm)");
current_axis = gca;
current_axis.YAxis.Exponent = 0;

subplot_handles(3) = subplot(3,1,3);
hold on
for quantile_index = 1:floor(number_of_quantiles/2)
    upper = saturation_state_quantiles(:,100-quantile_index);
    lower = saturation_state_quantiles(:,quantile_index);
    patch([interpolation_ages(age_location);flipud(interpolation_ages(age_location))],[lower;flipud(upper)],all_pH.colour,'FaceAlpha',0.026,'EdgeColor','None');
end
plot(interpolation_ages(age_location),saturation_state_quantiles(:,ceil(number_of_quantiles/2)),'-','Color',all_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);

ylim([0,12]);
set(gca,'YTick',0:2:12,'YTickLabels',num2str((0:2:12)'));


set(gca,'XDir','Reverse');
xlabel("Age (Ma)");
ylabel("Saturation State");


linkaxes(subplot_handles,'x');
xlim([min(interpolation_ages),max(interpolation_ages)]);

set(subplot_handles(1:2),'XTick',[]);
set(subplot_handles,'XDir','Reverse');

current_position = get(gcf,'Position');
set(gcf,'Position',[current_position(1:2),600,650])

exportgraphics(gcf,"./../../../Figures/d13C_pH_CO2_SaturationState_Evolutions.png","Resolution",600);
