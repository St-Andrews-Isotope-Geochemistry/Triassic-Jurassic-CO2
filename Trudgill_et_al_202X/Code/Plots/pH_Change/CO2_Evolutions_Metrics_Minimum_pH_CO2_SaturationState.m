clear

%% Load data
data_directory = "./../../../Data/";
boron_data = readtable(data_directory+"/Boron/TJ_d11B.xlsx","Sheet","Delta_Temperature");

posterior = jsondecode(fileread(data_directory+"/pH_change/posterior.json")).posterior;

interpolation_ages = jsondecode(fileread(data_directory+"/pH_change/posterior.json")).age;

raw_evolutions_metrics = readmatrix(data_directory+"/pH_Change/TJ_CO2_Evolutions_Metrics.csv");

pH_quantiles = raw_evolutions_metrics(1:3,:)';
saturation_state_quantiles = raw_evolutions_metrics(4:6,:)';
co2_quantiles = raw_evolutions_metrics(7:9,:)';

high_initial_pH.pH_quantiles = raw_evolutions_metrics(10:12,:)';
high_initial_pH.saturation_state_quantiles = raw_evolutions_metrics(13:15,:)';
high_initial_pH.co2_quantiles = raw_evolutions_metrics(16:18,:)';

raw_minimum_metrics = jsondecode(fileread(data_directory+"/Minimum_pH_Change/Metrics.json"));
raw_evolutions_metrics_delta = readmatrix(data_directory+"/pH_Change/Metrics_Delta.csv");

delta_pH_quantiles = raw_evolutions_metrics_delta(1:3,:)';
delta_saturation_state_quantiles = raw_evolutions_metrics_delta(4:6,:)';
delta_co2_quantiles = raw_evolutions_metrics_delta(7:9,:)';

delta_high_initial_pH.pH_quantiles = raw_evolutions_metrics_delta(10:12,:)';
delta_high_initial_pH.saturation_state_quantiles = raw_evolutions_metrics_delta(13:15,:)';
delta_high_initial_pH.co2_quantiles = raw_evolutions_metrics_delta(16:18,:)';

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

subplot_handles(1) = subplot(3,1,1);
hold on
patch([interpolation_ages;flipud(interpolation_ages)],[high_initial_pH.pH_quantiles(:,3);flipud(high_initial_pH.pH_quantiles(:,1))],high_initial_pH.colour,'FaceAlpha',0.2,'EdgeColor','None');
patch([interpolation_ages;flipud(interpolation_ages)],[pH_quantiles(:,3);flipud(pH_quantiles(:,1))],all_pH.colour,'FaceAlpha',0.3,'EdgeColor','None');

pH_handles(2) = plot(interpolation_ages,high_initial_pH.pH_quantiles(:,2),'-','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width/1.5);
pH_handles(1) = plot(interpolation_ages,pH_quantiles(:,2),'-','Color',all_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);


patch([initial_age,fliplr(initial_age)],[raw_minimum_metrics.initial.pH_025,raw_minimum_metrics.initial.pH_025,raw_minimum_metrics.initial.pH_975,raw_minimum_metrics.initial.pH_975],[0.7,0.7,0.7],'FaceAlpha',0.8,'EdgeColor','None');
plot(initial_age,[raw_minimum_metrics.initial.pH_median,raw_minimum_metrics.initial.pH_median],'Color','k','LineWidth',1.5);

patch([perturbation_age,fliplr(perturbation_age)],[raw_minimum_metrics.after.pH_025,raw_minimum_metrics.after.pH_025,raw_minimum_metrics.after.pH_975,raw_minimum_metrics.after.pH_975],[0.7,0.7,0.7],'FaceAlpha',0.8,'EdgeColor','None');
plot(perturbation_age,[raw_minimum_metrics.after.pH_median,raw_minimum_metrics.after.pH_median],'Color','k','LineWidth',2);

set(gca,'XDir','Reverse','XTick',[]);
ylabel("pH");

legend_handle = legend(pH_handles,["All evolutions","pH_i > "+num2str(round(raw_minimum_metrics.initial.pH_median-0.1,3))],'Location','SouthWest');

subplot_handles(2) = subplot(3,1,2);
hold on
patch([interpolation_ages;flipud(interpolation_ages)],[(delta_co2_quantiles(:,3));flipud((delta_co2_quantiles(:,1)))],all_pH.colour,'FaceAlpha',0.3,'EdgeColor','None');
patch([interpolation_ages;flipud(interpolation_ages)],[(delta_high_initial_pH.co2_quantiles(:,3));flipud((delta_high_initial_pH.co2_quantiles(:,1)))],high_initial_pH.colour,'FaceAlpha',0.2,'EdgeColor','None');

plot(interpolation_ages,(delta_high_initial_pH.co2_quantiles(:,2)),'-','Color',high_initial_pH.colour,'LineWidth',line_width/1.5);
plot(interpolation_ages,(delta_co2_quantiles(:,2)),'-','Color',all_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);

set(gca,'XTick',[]);
ylabel("\DeltaCO_2 (doublings)");

subplot_handles(3) = subplot(3,1,3);
hold on
patch([interpolation_ages;flipud(interpolation_ages)],[saturation_state_quantiles(:,3);flipud(saturation_state_quantiles(:,1))],all_pH.colour,'FaceAlpha',0.3,'EdgeColor','None');
patch([interpolation_ages;flipud(interpolation_ages)],[high_initial_pH.saturation_state_quantiles(:,3);flipud(high_initial_pH.saturation_state_quantiles(:,1))],high_initial_pH.colour,'FaceAlpha',0.2,'EdgeColor','None');

plot(interpolation_ages,high_initial_pH.saturation_state_quantiles(:,2),'-','Color',high_initial_pH.colour,'LineWidth',line_width/1.5);
plot(interpolation_ages,saturation_state_quantiles(:,2),'-','Color',all_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);

ylim([0,12]);
set(gca,'XDir','Reverse','YTick',[0:2:12],'YTickLabels',num2str([0:2:12]'));


set(gca,'XDir','Reverse');
xlabel("Age (Ma)");
ylabel("Saturation State");


linkaxes(subplot_handles,'x');
xlim([min(interpolation_ages),max(interpolation_ages)]);
set(subplot_handles,'XDir','Reverse');

exportgraphics(gcf,"./../../../Figures/d13C_pH_CO2_SaturationState_Evolutions.png");
