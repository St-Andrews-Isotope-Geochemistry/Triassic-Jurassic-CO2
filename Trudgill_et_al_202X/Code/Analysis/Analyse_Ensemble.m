clear

%% Load data
d18O_d13C = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Matlab");
d18O_d13C_averaged = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Averaged");

boron_data = readtable("./../../Data/TJ_d11B_pH.xlsx");
boron_data.age = boron_data.absolute_age;

raw_evolutions = readmatrix("./../../Data/TJ_CO2_Evolutions.csv");
reshaped_evolutions = reshape(raw_evolutions,[22,11,100000]);

evolutions.pH = squeeze(reshaped_evolutions(:,1,:));
evolutions.co2 = squeeze(reshaped_evolutions(:,2,:));
evolutions.saturation_state = squeeze(reshaped_evolutions(:,3,:));
evolutions.dic = squeeze(reshaped_evolutions(:,4,:));
evolutions.alkalinity = squeeze(reshaped_evolutions(:,5,:));
evolutions.temperature = squeeze(reshaped_evolutions(:,6,:));
evolutions.d11B = squeeze(reshaped_evolutions(:,7,:));
evolutions.calcium = squeeze(reshaped_evolutions(:,8,:));
evolutions.magnesium = squeeze(reshaped_evolutions(:,9,:));
evolutions.epsilon = squeeze(reshaped_evolutions(:,10,:));
evolutions.d11B_sw = squeeze(reshaped_evolutions(:,1,:));

clear raw_evolutions reshaped_evolutions

%% Get the initial subsample
evolutions.subsample_boolean = repmat(evolutions.saturation_state(1,:)>=5 & evolutions.saturation_state(1,:)<=10.7 & evolutions.co2(1,:)>=400 & evolutions.co2(1,:)<=5000 & all(evolutions.co2<10000) & all(evolutions.co2>0) & all(evolutions.saturation_state<12),size(evolutions.pH,1),1);
evolutions.pH_subsample = reshape(evolutions.pH(evolutions.subsample_boolean),22,[]);
evolutions.co2_subsample =  reshape(evolutions.co2(evolutions.subsample_boolean),22,[]);
evolutions.saturation_state_subsample = reshape(evolutions.saturation_state(evolutions.subsample_boolean),22,[]);

%% High pH subsample
high_initial_co2.boolean = repmat(evolutions.co2_subsample(1,:)>=2000,size(evolutions.pH_subsample,1),1);
high_initial_co2.pH = reshape(evolutions.pH_subsample(high_initial_co2.boolean),22,[]);
high_initial_co2.co2 =  reshape(evolutions.co2_subsample(high_initial_co2.boolean),22,[]);
high_initial_co2.saturation_state = reshape(evolutions.saturation_state_subsample(high_initial_co2.boolean),22,[]);

for distribution_index = 1:size(high_initial_co2.pH,1)
    high_initial_co2.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,high_initial_co2.pH(distribution_index,:));
    high_initial_co2.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,high_initial_co2.co2(distribution_index,:));
    high_initial_co2.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,high_initial_co2.saturation_state(distribution_index,:));
end

high_initial_co2.pH_quantiles = [high_initial_co2.pH_distributions.quantile(0.025),high_initial_co2.pH_distributions.quantile(0.5),high_initial_co2.pH_distributions.quantile(0.975)];
high_initial_co2.saturation_state_quantiles = [high_initial_co2.saturation_state_distributions.quantile(0.025),high_initial_co2.saturation_state_distributions.quantile(0.5),high_initial_co2.saturation_state_distributions.quantile(0.975)];
high_initial_co2.co2_quantiles = [high_initial_co2.co2_distributions.quantile(0.025),high_initial_co2.co2_distributions.quantile(0.5),high_initial_co2.co2_distributions.quantile(0.975)];

%% Medium CO2 subsample
medium_initial_co2.boolean = repmat(evolutions.co2_subsample(1,:)>=600 & evolutions.co2_subsample(1,:)<=2000,size(evolutions.pH_subsample,1),1);
medium_initial_co2.pH = reshape(evolutions.pH_subsample(medium_initial_co2.boolean),22,[]);
medium_initial_co2.co2 =  reshape(evolutions.co2_subsample(medium_initial_co2.boolean),22,[]);
medium_initial_co2.saturation_state = reshape(evolutions.saturation_state_subsample(medium_initial_co2.boolean),22,[]);

for distribution_index = 1:size(medium_initial_co2.pH,1)
    medium_initial_co2.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,medium_initial_co2.pH(distribution_index,:));
    medium_initial_co2.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,medium_initial_co2.co2(distribution_index,:));
    medium_initial_co2.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,medium_initial_co2.saturation_state(distribution_index,:));
end

medium_initial_co2.pH_quantiles = [medium_initial_co2.pH_distributions.quantile(0.025),medium_initial_co2.pH_distributions.quantile(0.5),medium_initial_co2.pH_distributions.quantile(0.975)];
medium_initial_co2.saturation_state_quantiles = [medium_initial_co2.saturation_state_distributions.quantile(0.025),medium_initial_co2.saturation_state_distributions.quantile(0.5),medium_initial_co2.saturation_state_distributions.quantile(0.975)];
medium_initial_co2.co2_quantiles = [medium_initial_co2.co2_distributions.quantile(0.025),medium_initial_co2.co2_distributions.quantile(0.5),medium_initial_co2.co2_distributions.quantile(0.975)];

%% Low pH subsample
low_initial_co2.boolean = repmat(evolutions.co2_subsample(1,:)>400 & evolutions.co2_subsample(1,:)<600,size(evolutions.pH_subsample,1),1);
low_initial_co2.pH = reshape(evolutions.pH_subsample(low_initial_co2.boolean),22,[]);
low_initial_co2.co2 =  reshape(evolutions.co2_subsample(low_initial_co2.boolean),22,[]);
low_initial_co2.saturation_state = reshape(evolutions.saturation_state_subsample(low_initial_co2.boolean),22,[]);

for distribution_index = 1:size(low_initial_co2.pH,1)
    low_initial_co2.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,low_initial_co2.pH(distribution_index,:));
    low_initial_co2.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,low_initial_co2.co2(distribution_index,:));
    low_initial_co2.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,low_initial_co2.saturation_state(distribution_index,:));
end

low_initial_co2.pH_quantiles = [low_initial_co2.pH_distributions.quantile(0.025),low_initial_co2.pH_distributions.quantile(0.5),low_initial_co2.pH_distributions.quantile(0.975)];
low_initial_co2.saturation_state_quantiles = [low_initial_co2.saturation_state_distributions.quantile(0.025),low_initial_co2.saturation_state_distributions.quantile(0.5),low_initial_co2.saturation_state_distributions.quantile(0.975)];
low_initial_co2.co2_quantiles = [low_initial_co2.co2_distributions.quantile(0.025),low_initial_co2.co2_distributions.quantile(0.5),low_initial_co2.co2_distributions.quantile(0.975)];

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

%% Extra plot for the distribution of CO2 change
co2_change = max(evolutions.co2_subsample(10:end,:))-mean(evolutions.co2_subsample (1:9,:));
co2_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(-1000:100:10000,co2_change).normalise();

co2_change_quantiles = [co2_change_distribution.quantile(0.05),co2_change_distribution.quantile(0.5),co2_change_distribution.quantile(0.95)];

figure(2);
clf
hold on
co2_change_distribution.plot();
plot([co2_change_quantiles;co2_change_quantiles],[0;1],'k');

ylim([0,max(co2_change_distribution.probabilities)]);
xlabel("CO_2 (ppm)");
ylabel("Probability");

%% pH
clf
hold on
high_initial_co2.pH_change = min(high_initial_co2.pH(10:end,:))-mean(high_initial_co2.pH(1:9,:));
high_initial_co2.pH_change_distritbuion = Geochemistry_Helpers.Distribution.fromSamples(-2:0.01:1,high_initial_co2.pH_change);
high_initial_co2.pH_change_distritbuion.plot();

high_initial_co2.pH_change_distritbuion.quantile(0.95);
high_initial_co2.pH_change_distritbuion.quantile(0.5);
high_initial_co2.pH_change_distritbuion.standard_deviation();

low_initial_co2.pH_change = min(low_initial_co2.pH(10:end,:))-mean(low_initial_co2.pH(1:9,:));
low_initial_co2.pH_change_distritbuion = Geochemistry_Helpers.Distribution.fromSamples(-2:0.01:1,low_initial_co2.pH_change);
low_initial_co2.pH_change_distritbuion.plot();

low_initial_co2.pH_change_distritbuion.quantile(0.5);
low_initial_co2.pH_change_distritbuion.standard_deviation();

pH_changes = min(evolutions.pH_subsample(10:end,:))-mean(evolutions.pH_subsample(1:9,:));
pH_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(-2:0.01:1,pH_changes);
pH_change_distribution.plot();

pH_change_distribution.quantile(0.95);
pH_change_distribution.quantile(0.5);
pH_change_distribution.standard_deviation()*2