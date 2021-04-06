clear

%% Load data
boron_data = readtable("./../Data/TJ_d11B_pH.xlsx");
boron_data.age = boron_data.absolute_age;

raw_evolutions = readmatrix("./../Data/TJ_CO2_Evolutions.csv");
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
high_initial_pH.boolean = repmat(evolutions.pH_subsample(1,:)>=8,size(evolutions.pH_subsample,1),1);
high_initial_pH.pH = reshape(evolutions.pH_subsample(high_initial_pH.boolean),22,[]);
high_initial_pH.co2 =  reshape(evolutions.co2_subsample(high_initial_pH.boolean),22,[]);
high_initial_pH.saturation_state = reshape(evolutions.saturation_state_subsample(high_initial_pH.boolean),22,[]);

for distribution_index = 1:size(high_initial_pH.pH,1)
    high_initial_pH.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,high_initial_pH.pH(distribution_index,:));
    high_initial_pH.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,high_initial_pH.co2(distribution_index,:));
    high_initial_pH.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,high_initial_pH.saturation_state(distribution_index,:));
end

high_initial_pH.pH_quantiles = [high_initial_pH.pH_distributions.quantile(0.025),high_initial_pH.pH_distributions.quantile(0.5),high_initial_pH.pH_distributions.quantile(0.975)];
high_initial_pH.saturation_state_quantiles = [high_initial_pH.saturation_state_distributions.quantile(0.025),high_initial_pH.saturation_state_distributions.quantile(0.5),high_initial_pH.saturation_state_distributions.quantile(0.975)];
high_initial_pH.co2_quantiles = [high_initial_pH.co2_distributions.quantile(0.025),high_initial_pH.co2_distributions.quantile(0.5),high_initial_pH.co2_distributions.quantile(0.975)];

%% Low pH subsample
low_initial_pH.boolean = repmat(evolutions.pH_subsample(1,:)<8,size(evolutions.pH_subsample,1),1);
low_initial_pH.pH = reshape(evolutions.pH_subsample(low_initial_pH.boolean),22,[]);
low_initial_pH.co2 =  reshape(evolutions.co2_subsample(low_initial_pH.boolean),22,[]);
low_initial_pH.saturation_state = reshape(evolutions.saturation_state_subsample(low_initial_pH.boolean),22,[]);

for distribution_index = 1:size(low_initial_pH.pH,1)
    low_initial_pH.pH_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,low_initial_pH.pH(distribution_index,:));
    low_initial_pH.co2_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,low_initial_pH.co2(distribution_index,:));
    low_initial_pH.saturation_state_distributions(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples(0:0.1:20,low_initial_pH.saturation_state(distribution_index,:));
end

low_initial_pH.pH_quantiles = [low_initial_pH.pH_distributions.quantile(0.025),low_initial_pH.pH_distributions.quantile(0.5),low_initial_pH.pH_distributions.quantile(0.975)];
low_initial_pH.saturation_state_quantiles = [low_initial_pH.saturation_state_distributions.quantile(0.025),low_initial_pH.saturation_state_distributions.quantile(0.5),low_initial_pH.saturation_state_distributions.quantile(0.975)];
low_initial_pH.co2_quantiles = [low_initial_pH.co2_distributions.quantile(0.025),low_initial_pH.co2_distributions.quantile(0.5),low_initial_pH.co2_distributions.quantile(0.975)];

%% Make a plot like figure 2 in the paper
high_initial_pH.colour = [0,0,0.6];
low_initial_pH.colour = [0,0.6,0.2];
marker_size = 4;
line_width = 1.2;

figure(1);
clf
subplot(3,1,1);
hold on
high_initial_pH.handles{1} = plot(boron_data.age,high_initial_pH.pH_quantiles(:,3),'^','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
high_initial_pH.handles{2} = plot(boron_data.age,high_initial_pH.pH_quantiles(:,2),'o-','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
high_initial_pH.handles{3} = plot(boron_data.age,high_initial_pH.pH_quantiles(:,1),'v','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);

low_initial_pH.handles{1} = plot(boron_data.age,low_initial_pH.pH_quantiles(:,3),'^','Color',low_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
low_initial_pH.handles{2} = plot(boron_data.age,low_initial_pH.pH_quantiles(:,2),'o-','Color',low_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
low_initial_pH.handles{3} = plot(boron_data.age,low_initial_pH.pH_quantiles(:,1),'v','Color',low_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);

set(gca,'XDir','Reverse','XTick',[]);
ylabel("pH");

legend_handle_1 = legend([high_initial_pH.handles{:}],["Upper 95%","Median","Lower 95%"],'Location','SouthWest');

dummy_axis = axes('position',get(gca,'position'),'visible','off');
legend_handle_2 = legend(dummy_axis,[low_initial_pH.handles{:}],["Upper 95%","Median","Lower 95%"],'Location','SouthWest');

title(legend_handle_2,"pH_i <8.0");
title(legend_handle_1,"pH_i >8.0");

legend_handle_1.Position = legend_handle_2.Position+[legend_handle_2.Position(3)+0.01,0,0,0];

subplot(3,1,2);
hold on
plot(boron_data.age,high_initial_pH.saturation_state_quantiles(:,3),'^','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
plot(boron_data.age,high_initial_pH.saturation_state_quantiles(:,2),'o-','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
plot(boron_data.age,high_initial_pH.saturation_state_quantiles(:,1),'v','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);

plot(boron_data.age,low_initial_pH.saturation_state_quantiles(:,3),'^','Color',low_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
plot(boron_data.age,low_initial_pH.saturation_state_quantiles(:,2),'o-','Color',low_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
plot(boron_data.age,low_initial_pH.saturation_state_quantiles(:,1),'v','Color',low_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);

set(gca,'XDir','Reverse','XTick',[]);
ylabel("Saturation State");

subplot(3,1,3);
hold on
plot(boron_data.age,log2(high_initial_pH.co2_quantiles(:,3)),'^','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
plot(boron_data.age,log2(high_initial_pH.co2_quantiles(:,2)),'o-','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
plot(boron_data.age,log2(high_initial_pH.co2_quantiles(:,1)),'v','Color',high_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);

plot(boron_data.age,log2(low_initial_pH.co2_quantiles(:,3)),'^','Color',low_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
plot(boron_data.age,log2(low_initial_pH.co2_quantiles(:,2)),'o-','Color',low_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);
plot(boron_data.age,log2(low_initial_pH.co2_quantiles(:,1)),'v','Color',low_initial_pH.colour,'MarkerSize',marker_size,'LineWidth',line_width);

set(gca,'XDir','Reverse');
xlabel("Age (Ma)");
ylabel("CO_2 (ppm)");

current_axis = gca;
current_axis.YAxis.MinorTickValues = log2([400:400:12800*2]);
y_tick_locations = 400*2.^[-1,0,1,2,3,4,5];
set(gca,'YTick',log2(y_tick_locations),'YTickLabels',num2str(y_tick_locations'),'YMinorTick','On');

%% Extra plot for the distribution of CO2 change
co2_change = max(evolutions.co2_subsample (10:end,:))-mean(evolutions.co2_subsample (1:9,:));
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
