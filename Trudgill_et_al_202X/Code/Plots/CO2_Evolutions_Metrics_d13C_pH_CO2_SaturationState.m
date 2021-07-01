clear

%% Load data
d18O_d13C = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab");
d18O_d13C_averaged = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Averaged");

boron_data = readtable("./../../Data/TJ_d11B_pH.xlsx");
boron_data.age = boron_data.absolute_age;

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
d18O_d13C.d13C(d18O_d13C.absolute_age<min(boron_data.age)) = NaN;
d18O_d13C_averaged.d13C(d18O_d13C_averaged.absolute_age<min(boron_data.age)) = NaN;

%% Make a plot like figure 2 in the paper
high_initial_co2.colour = [0,0,0.6];
medium_initial_co2.colour = [0.6,0,0.2];
low_initial_co2.colour = [0,0.6,0.2];
marker_size = 4;
line_width = 1.2;

figure(1);
clf
subplot_handles(1) = subplot(4,1,1);
hold on
plot(d18O_d13C.absolute_age,d18O_d13C.d13C,'ok','MarkerFaceColor','k','MarkerSize',marker_size/2,'LineWidth',line_width);
% plot(d18O_d13C_averaged.absolute_age,d18O_d13C_averaged.d13C,'k-','LineWidth',line_width);

ylim([1,5]);
set(gca,'XDir','Reverse','XTick',[]);
ylabel("\delta^{13}C (‰)");

box off;

subplot_handles(2) = subplot(4,1,2);
hold on
% high_initial_co2.handles{1} = plot(boron_data.age,high_initial_co2.pH_quantiles(:,3),'^','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);
% high_initial_co2.handles{2} = plot(boron_data.age,high_initial_co2.pH_quantiles(:,2),'o-','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);
% high_initial_co2.handles{3} = plot(boron_data.age,high_initial_co2.pH_quantiles(:,1),'v','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);

% medium_initial_co2.handles{1} = plot(boron_data.age,medium_initial_co2.pH_quantiles(:,3),'^','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);
medium_initial_co2.handles{1} = plot(boron_data.age,medium_initial_co2.pH_quantiles(:,2),'o-','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
medium_initial_co2.handles{2} = plot(boron_data.age,medium_initial_co2.pH_quantiles(:,1),'v','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width/1.5);

low_initial_co2.handles{1} = plot(boron_data.age,low_initial_co2.pH_quantiles(:,3),'^','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width/1.5);
low_initial_co2.handles{2} = plot(boron_data.age,low_initial_co2.pH_quantiles(:,2),'o-','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',low_initial_co2.colour);
% low_initial_co2.handles{3} = plot(boron_data.age,low_initial_co2.pH_quantiles(:,1),'v','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);


set(gca,'XDir','Reverse','XTick',[]);
ylabel("pH");

% legend_handle_1 = legend([medium_initial_co2.handles{:}],["Upper 95%","Median","Lower 95%"],'Location','SouthWest');
legend_handle_1 = legend([low_initial_co2.handles{2},medium_initial_co2.handles{1}],["CO_2_i = 400-600ppm","CO_2_i = 600-2000ppm"],'Location','SouthWest');

% dummy_axis = axes('position',get(gca,'position'),'visible','off');
% legend_handle_2 = legend(dummy_axis,[low_initial_co2.handles{:}],["Upper 95%","Median","Lower 95%"],'Location','SouthWest');

% title(legend_handle_2,"CO_2_i 400-600ppm");
% title(legend_handle_1,"CO_2_i 600-2000ppm");

% legend_handle_1.Position = legend_handle_2.Position+[legend_handle_2.Position(3)+0.01,0,0,0];

subplot_handles(3) = subplot(4,1,3);
hold on
% plot(boron_data.age,log2(high_initial_co2.co2_quantiles(:,3)),'^','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);
% plot(boron_data.age,log2(high_initial_co2.co2_quantiles(:,2)),'o-','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);
% plot(boron_data.age,log2(high_initial_co2.co2_quantiles(:,1)),'v','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);

plot(boron_data.age,log2(medium_initial_co2.co2_quantiles(:,3)),'^','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width/1.5);
plot(boron_data.age,log2(medium_initial_co2.co2_quantiles(:,2)),'o-','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
% plot(boron_data.age,log2(medium_initial_co2.co2_quantiles(:,1)),'v','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);

% plot(boron_data.age,log2(low_initial_co2.co2_quantiles(:,3)),'^','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);
plot(boron_data.age,log2(low_initial_co2.co2_quantiles(:,2)),'o-','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',low_initial_co2.colour);
plot(boron_data.age,log2(low_initial_co2.co2_quantiles(:,1)),'v','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width/1.5);

set(gca,'XDir','Reverse','XTick',[]);
ylabel("CO_2 (ppm)");

current_axis = gca;
current_axis.YAxis.MinorTickValues = log2([400:400:12800*2]);
y_tick_locations = 400*2.^[-1,0,1,2,3,4,5];
set(gca,'YTick',log2(y_tick_locations),'YTickLabels',num2str(y_tick_locations'),'YMinorTick','On');

subplot_handles(4) = subplot(4,1,4);
hold on
% plot(boron_data.age,high_initial_co2.saturation_state_quantiles(:,3),'^','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);
% plot(boron_data.age,high_initial_co2.saturation_state_quantiles(:,2),'o-','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);
% plot(boron_data.age,high_initial_co2.saturation_state_quantiles(:,1),'v','Color',high_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width);

plot(boron_data.age,medium_initial_co2.saturation_state_quantiles(:,3),'^','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width/1.5);
plot(boron_data.age,medium_initial_co2.saturation_state_quantiles(:,2),'o-','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',medium_initial_co2.colour);
plot(boron_data.age,medium_initial_co2.saturation_state_quantiles(:,1),'v','Color',medium_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width/1.5);

plot(boron_data.age,low_initial_co2.saturation_state_quantiles(:,3),'^','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width/1.5);
plot(boron_data.age,low_initial_co2.saturation_state_quantiles(:,2),'o-','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width,'MarkerFaceColor',low_initial_co2.colour);
plot(boron_data.age,low_initial_co2.saturation_state_quantiles(:,1),'v','Color',low_initial_co2.colour,'MarkerSize',marker_size,'LineWidth',line_width/1.5);

set(gca,'XDir','Reverse');
xlabel("Age (Ma)");
ylabel("Saturation State");

original_xlim = xlim;

linkaxes(subplot_handles,'x');
xlim(original_xlim);

exportgraphics(gcf,"./../../Figures/d13C_pH_CO2_SaturationState_Evolutions.png");
