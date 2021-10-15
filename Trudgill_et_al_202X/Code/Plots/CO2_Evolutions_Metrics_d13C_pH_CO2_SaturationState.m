clear

%% Load data
% d18O_d13C = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab");
% d18O_d13C_averaged = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Averaged");

boron_data = readtable("./../../Data/TJ_d11B.xlsx","Sheet","Delta_Temperature");
% interpolation_ages = boron_data.absolute_age;
interpolation_ages = unique(sort([boron_data.age',linspace(min(boron_data.age),max(boron_data.age),80)]));

raw_evolutions_metrics = readmatrix("./../../Data/TJ_CO2_Evolutions_Metrics.csv");

low_initial_co2.pH_quantiles = raw_evolutions_metrics(1:3,:)';
medium_initial_co2.pH_quantiles = raw_evolutions_metrics(4:6,:)';
high_initial_co2.pH_quantiles = raw_evolutions_metrics(7:9,:)';

low_initial_co2.saturation_state_quantiles = raw_evolutions_metrics(10:12,:)';
medium_initial_co2.saturation_state_quantiles = raw_evolutions_metrics(13:15,:)';
high_initial_co2.saturation_state_quantiles = raw_evolutions_metrics(16:18,:)';

low_initial_co2.co2_quantiles = raw_evolutions_metrics(19:21,:)';
medium_initial_co2.co2_quantiles = raw_evolutions_metrics(22:24,:)';
high_initial_co2.co2_quantiles = raw_evolutions_metrics(25:27,:)';

%%
% d18O_d13C.d13C(d18O_d13C.absolute_age<min(interpolation_ages)) = NaN;
% d18O_d13C_averaged.d13C(d18O_d13C_averaged.absolute_age<min(interpolation_ages)) = NaN;

%% Make a plot like figure 2 in the paper
high_initial_co2.colour = [0.6,0,0];
medium_initial_co2.colour = [0.8,0.8,0.2];
low_initial_co2.colour = [0,0,0.6];
marker_size = 4;
line_width = 1.2;

figure(1);
clf

subplot_handles(1) = subplot(3,1,1);
hold on
patch([interpolation_ages,fliplr(interpolation_ages)],[high_initial_co2.pH_quantiles(:,3);flipud(high_initial_co2.pH_quantiles(:,1))],high_initial_co2.colour,'FaceAlpha',0.2,'EdgeColor','None');
% patch([interpolation_ages,fliplr(interpolation_ages)],[log2(medium_initial_co2.pH_quantiles(:,3));flipud(log2(medium_initial_co2.pH_quantiles(:,1)))],medium_initial_co2.colour,'FaceAlpha',0.2,'EdgeColor','None');
patch([interpolation_ages,fliplr(interpolation_ages)],[low_initial_co2.pH_quantiles(:,3);flipud(low_initial_co2.pH_quantiles(:,1))],low_initial_co2.colour,'FaceAlpha',0.2,'EdgeColor','None');

high_initial_co2.handles{1} = plot(interpolation_ages,high_initial_co2.pH_quantiles(:,2),'-','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
% plot(interpolation_ages,log2(medium_initial_co2.pH_quantiles(:,2)),'-','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
low_initial_co2.handles{1} = plot(interpolation_ages,low_initial_co2.pH_quantiles(:,2),'-','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',low_initial_co2.colour);

set(gca,'XDir','Reverse','XTick',[]);
ylabel("pH");

% legend_handle_1 = legend([medium_initial_co2.handles{:}],["Upper 95%","Median","Lower 95%"],'Location','SouthWest');
legend_handle_1 = legend([low_initial_co2.handles{1},high_initial_co2.handles{1}],["CO_2_i < 800ppm","CO_2_i >= 800ppm"],'Location','SouthWest');

% dummy_axis = axes('position',get(gca,'position'),'visible','off');
% legend_handle_2 = legend(dummy_axis,[low_initial_co2.handles{:}],["Upper 95%","Median","Lower 95%"],'Location','SouthWest');

% title(legend_handle_2,"CO_2_i 400-600ppm");
% title(legend_handle_1,"CO_2_i 600-2000ppm");

% legend_handle_1.Position = legend_handle_2.Position+[legend_handle_2.Position(3)+0.01,0,0,0];

subplot_handles(2) = subplot(3,1,2);
hold on
patch([interpolation_ages,fliplr(interpolation_ages)],[log2(high_initial_co2.co2_quantiles(:,3));flipud(log2(high_initial_co2.co2_quantiles(:,1)))],high_initial_co2.colour,'FaceAlpha',0.2,'EdgeColor','None');
% patch([interpolation_ages,fliplr(interpolation_ages)],[log2(medium_initial_co2.co2_quantiles(:,3));flipud(log2(medium_initial_co2.co2_quantiles(:,1)))],medium_initial_co2.colour,'FaceAlpha',0.2,'EdgeColor','None');
patch([interpolation_ages,fliplr(interpolation_ages)],[log2(low_initial_co2.co2_quantiles(:,3));flipud(log2(low_initial_co2.co2_quantiles(:,1)))],low_initial_co2.colour,'FaceAlpha',0.2,'EdgeColor','None');

plot(interpolation_ages,log2(high_initial_co2.co2_quantiles(:,2)),'-','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
% plot(interpolation_ages,log2(medium_initial_co2.co2_quantiles(:,2)),'-','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
plot(interpolation_ages,log2(low_initial_co2.co2_quantiles(:,2)),'-','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',low_initial_co2.colour);

set(gca,'XDir','Reverse','XTick',[]);
ylabel("CO_2 (ppm)");

current_axis = gca;
current_axis.YAxis.MinorTickValues = log2(400*2.^[-1,0,1,2,3,4,5,6,7,8]);
y_tick_locations = 400*2.^[-1,0,1,2,3,4,5,6,7,8];
set(gca,'YTick',log2(y_tick_locations),'YTickLabels',num2str(y_tick_locations'),'YMinorTick','On');

subplot_handles(3) = subplot(3,1,3);
hold on
patch([interpolation_ages,fliplr(interpolation_ages)],[high_initial_co2.saturation_state_quantiles(:,3);flipud(high_initial_co2.saturation_state_quantiles(:,1))],high_initial_co2.colour,'FaceAlpha',0.2,'EdgeColor','None');
% patch([interpolation_ages,fliplr(interpolation_ages)],[log2(medium_initial_co2.saturation_state_quantiles(:,3));flipud(log2(medium_initial_co2.saturation_state_quantiles(:,1)))],medium_initial_co2.colour,'FaceAlpha',0.2,'EdgeColor','None');
patch([interpolation_ages,fliplr(interpolation_ages)],[low_initial_co2.saturation_state_quantiles(:,3);flipud(low_initial_co2.saturation_state_quantiles(:,1))],low_initial_co2.colour,'FaceAlpha',0.2,'EdgeColor','None');

plot(interpolation_ages,high_initial_co2.saturation_state_quantiles(:,2),'-','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
% plot(interpolation_ages,log2(medium_initial_co2.saturation_state_quantiles(:,2)),'-','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
plot(interpolation_ages,low_initial_co2.saturation_state_quantiles(:,2),'-','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',low_initial_co2.colour);

ylim([0,12]);
set(gca,'XDir','Reverse','YTick',[0:2:12],'YTickLabels',num2str([0:2:12]'));


set(gca,'XDir','Reverse');
xlabel("Age (Ma)");
ylabel("Saturation State");

original_xlim = xlim;

linkaxes(subplot_handles,'x');
xlim(original_xlim);

exportgraphics(gcf,"./../../Figures/d13C_pH_CO2_SaturationState_Evolutions.png");
