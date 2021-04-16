clear

%% Load data
data = jsondecode(fileread("./../../Data/Minimum_pH_Variation.json"));

independent = containers.Map();
independent("calcium") = data{1}.calcium;
independent("magnesium") = data{2}.magnesium;
independent("temperature_change") = data{3}.temperature_change;
independent("initial_temperature") = data{4}.initial_temperature;
independent("saturation_state") = data{5}.saturation_state;
independent("co2") = data{6}.co2;
independent("epsilon") = data{7}.epsilon;

pH_difference = containers.Map();
pH_difference("calcium") = data{1}.pH_difference;
pH_difference("magnesium") = data{2}.pH_difference;
pH_difference("temperature_change") = data{3}.pH_difference;
pH_difference("initial_temperature") = data{4}.pH_difference;
pH_difference("saturation_state") = data{5}.pH_difference;
pH_difference("co2") = data{6}.pH_difference;
pH_difference("epsilon") = data{7}.pH_difference;

colour_map = get(gca,'ColorOrder');
% colour_map = jet(6);
colours = containers.Map();
colours("calcium") = colour_map(1,:);
colours("magnesium") = colour_map(2,:);
colours("temperature_change") = colour_map(3,:);
colours("initial_temperature") = colour_map(4,:);
colours("saturation_state") = colour_map(5,:);
colours("co2") = colour_map(6,:);
colours("epsilon") = colour_map(7,:);

%%
pH_difference_limits = [-0.8,-0.2];

clf
ca_axis = axes();
plot_handle = plot(independent("calcium"),pH_difference("calcium"),'Color',colours("calcium"),'LineWidth',2,'Parent',ca_axis);
xlim([min(independent("calcium")),max(independent("calcium"))]);

ylabel("\DeltapH");

mg_axis = axes();
plot(independent("magnesium"),pH_difference("magnesium"),'Color',colours("magnesium"),'LineWidth',2,'Parent',mg_axis);
xlim([min(independent("magnesium")),max(independent("magnesium"))]);

temperature_change_axis = axes();
plot(independent("temperature_change"),pH_difference("temperature_change"),'Color',colours("temperature_change"),'LineWidth',2,'Parent',temperature_change_axis);
xlim([min(independent("temperature_change")),max(independent("temperature_change"))]);

initial_temperature_axis = axes();
plot(independent("initial_temperature"),pH_difference("initial_temperature"),'Color',colours("initial_temperature"),'LineWidth',2,'Parent',initial_temperature_axis);
xlim([min(independent("initial_temperature")),max(independent("initial_temperature"))]);

saturation_state_axis = axes();
plot(independent("saturation_state"),pH_difference("saturation_state"),'Color',colours("saturation_state"),'LineWidth',2,'Parent',saturation_state_axis);
xlim([min(independent("saturation_state")),max(independent("saturation_state"))]);

co2_axis = axes();
plot(independent("co2")*1e6,pH_difference("co2"),'Color',colours("co2"),'LineWidth',2,'Parent',co2_axis);
xlim([min(independent("co2")),max(independent("co2"))]);

epsilon_axis = axes();
plot(independent("epsilon"),pH_difference("epsilon"),'Color',colours("epsilon"),'LineWidth',2,'Parent',epsilon_axis);
xlim([min(independent("epsilon")),max(independent("epsilon"))]);

set(ca_axis,'Position',ca_axis.Position+[0,0.04,0,0]);
original_position = get(ca_axis,"Position");
all_axes = [ca_axis,mg_axis,temperature_change_axis,initial_temperature_axis,saturation_state_axis,co2_axis,epsilon_axis];
for current_axis = 2:numel(all_axes)
    set(all_axes(current_axis),'Color','None','YTick',[]);
end
for current_axis = 1:numel(all_axes)
    set(all_axes(current_axis),'Position',original_position+[-0.03,0.08,-0.08,-0.16],'XTick',[]);
end
linkaxes(all_axes,'y');
ylim(pH_difference_limits);

new_position = get(ca_axis,'Position');

ca_xaxis = axes('Position',(new_position.*[1,1,1,1e-3])+[0,0,0,0],'Color','None');
plot(17,0,'o','Color',colours("calcium"),'MarkerFaceColor',colours("calcium"),'Parent',ca_xaxis);
set(ca_xaxis,'XColor',colours("calcium"));
current_label = xlabel("Calcium");
set(current_label,'Position',[max(independent("calcium"))+(max(independent("calcium"))-min(independent("calcium")))*0.02,12],'HorizontalAlignment','Left');
xlim([min(independent("calcium")),max(independent("calcium"))]);

mg_xaxis = axes('Position',(new_position.*[1,1,1,1e-3])+[0,-0.05,0,0],'Color','None');
plot(28,0,'o','Color',colours("magnesium"),'MarkerFaceColor',colours("magnesium"),'Parent',mg_xaxis);
set(mg_xaxis,'XColor',colours("magnesium"));
current_label = xlabel("Magnesium");
set(current_label,'Position',[max(independent("magnesium"))+(max(independent("magnesium"))-min(independent("magnesium")))*0.02,12],'HorizontalAlignment','Left');
xlim([min(independent("magnesium")),max(independent("magnesium"))]);

temperature_change_xaxis = axes('Position',(new_position.*[1,1,1,1e-3])+[0,-0.1,0,0],'Color','None');
plot(0,0,'o','Color',colours("temperature_change"),'MarkerFaceColor',colours("temperature_change"),'Parent',temperature_change_xaxis);
set(temperature_change_xaxis,'XColor',colours("temperature_change"));
current_label = xlabel("\DeltaTemperature",'Position',[21.8,12]);
set(current_label,'Position',[max(independent("temperature_change"))+(max(independent("temperature_change"))-min(independent("temperature_change")))*0.02,12],'HorizontalAlignment','Left');
xlim([min(independent("temperature_change")),max(independent("temperature_change"))]);

initial_temperature_xaxis =axes('Position',(new_position.*[1,1,1,1e-3])+[0,-0.15,0,0],'Color','None','XAxisLocation','Top');
plot(21.4,0,'o','Color',colours("initial_temperature"),'MarkerFaceColor',colours("initial_temperature"),'Parent',initial_temperature_xaxis);
set(initial_temperature_xaxis,'XColor',colours("initial_temperature"));
current_label = xlabel("Initial temperature");
set(current_label,'Position',[max(independent("initial_temperature"))+(max(independent("initial_temperature"))-min(independent("initial_temperature")))*0.02,12],'HorizontalAlignment','Left');
xlim([min(independent("initial_temperature")),max(independent("initial_temperature"))]);

epsilon_xaxis =axes('Position',(new_position.*[1,1,1,1e-3])+[0,new_position(4)+0,0,0],'Color','None','XAxisLocation','Top');
plot(27.2,0,'o','Color',colours("epsilon"),'MarkerFaceColor',colours("epsilon"),'Parent',epsilon_xaxis);
set(epsilon_xaxis,'XColor',colours("epsilon"));
current_label = xlabel("Epsilon");
set(current_label,'Position',[max(independent("epsilon"))+(max(independent("epsilon"))-min(independent("epsilon")))*0.02,20],'HorizontalAlignment','Left');
xlim([min(independent("epsilon")),max(independent("epsilon"))]);

saturation_state_xaxis = axes('Position',(new_position.*[1,1,1,1e-3])+[0,new_position(4)+0.05,0,0],'Color','None','XAxisLocation','Top');
plot(10.7,0,'o','Color',colours("saturation_state"),'MarkerFaceColor',colours("saturation_state"),'Parent',saturation_state_xaxis);
set(saturation_state_xaxis,'XColor',colours("saturation_state"));
current_label = xlabel("Saturation State");
set(current_label,'Position',[max(independent("saturation_state"))+(max(independent("saturation_state"))-min(independent("saturation_state")))*0.02,20],'HorizontalAlignment','Left');
xlim([min(independent("saturation_state")),max(independent("saturation_state"))]);

co2_xaxis = axes('Position',(new_position.*[1,1,1,1e-3])+[0,new_position(4)+0+0.1,0,0],'Color','None','XAxisLocation','Top');
plot(500,0,'o','Color',colours("co2"),'MarkerFaceColor',colours("co2"),'Parent',co2_xaxis);
set(co2_xaxis,'XColor',colours("co2"));
current_label = xlabel("CO_2");
set(current_label,'Position',[(max(independent("co2"))+(max(independent("co2"))-min(independent("co2")))*0.02)*1e6,20],'HorizontalAlignment','Left');
xlim([min(independent("co2")),max(independent("co2"))]*1e6);

exportgraphics(gcf,"./../../Figures/Delta_pH_Sensitivity.png");

