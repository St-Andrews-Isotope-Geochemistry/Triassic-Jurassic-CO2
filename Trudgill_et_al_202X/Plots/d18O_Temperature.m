clear
% data = readtable("./../Data/d18O_d13C","Sheet","Matlab");
averaged = readtable("./../Data/d18O_d13C","Sheet","Averaged");

%%
age_midpoints = averaged.absolute_age(1:end-1) + diff(averaged.absolute_age);
oneil_difference = diff(averaged.oneil_temperature);
kim_oneil_difference = diff(averaged.kim_oneil_temperature);
hansen_difference = diff(averaged.hansen_temperature);
anderson_arthur_difference = diff(averaged.anderson_arthur_temperature);

colours = get(gca,'colororder');
oneil_colour = colours(1,:);
kim_oneil_colour = colours(2,:);
hansen_colour = colours(3,:);
anderson_arthur_colour = colours(5,:);

%%
figure(1);
clf
hold on

for depth_index = 1:10
    plot([averaged.absolute_age(depth_index),averaged.absolute_age(depth_index)],averaged.oneil_temperature(depth_index)+[-2*averaged.oneil_temperature_uncertainty(depth_index),2*averaged.oneil_temperature_uncertainty(depth_index)],'Color',oneil_colour);
    plot([averaged.absolute_age(depth_index),averaged.absolute_age(depth_index)],averaged.kim_oneil_temperature(depth_index)+[-2*averaged.kim_oneil_temperature_uncertainty(depth_index),2*averaged.kim_oneil_temperature_uncertainty(depth_index)],'Color',kim_oneil_colour);
    plot([averaged.absolute_age(depth_index),averaged.absolute_age(depth_index)],averaged.hansen_temperature(depth_index)+[-2*averaged.hansen_temperature_uncertainty(depth_index),2*averaged.hansen_temperature_uncertainty(depth_index)],'Color',hansen_colour);
    plot([averaged.absolute_age(depth_index),averaged.absolute_age(depth_index)],averaged.anderson_arthur_temperature(depth_index)+[-2*averaged.anderson_arthur_temperature_uncertainty(depth_index),2*averaged.anderson_arthur_temperature_uncertainty(depth_index)],'Color',anderson_arthur_colour);
end
oneil_handle = plot(averaged.absolute_age(1:10),averaged.oneil_temperature(1:10),'.-','Color',oneil_colour);
kim_oneil_handle = plot(averaged.absolute_age(1:10),averaged.kim_oneil_temperature(1:10),'.-','Color',kim_oneil_colour);
hansen_handle = plot(averaged.absolute_age(1:10),averaged.hansen_temperature(1:10),'.-','Color',hansen_colour);
anderson_arthur_handle = plot(averaged.absolute_age(1:10),averaged.anderson_arthur_temperature(1:10),'.-','Color',anderson_arthur_colour);

for depth_index = 10:height(averaged)
    plot([averaged.absolute_age(depth_index),averaged.absolute_age(depth_index)],averaged.oneil_temperature(depth_index)+[-2*averaged.oneil_temperature_uncertainty(depth_index),2*averaged.oneil_temperature_uncertainty(depth_index)],'Color',oneil_colour);
    plot([averaged.absolute_age(depth_index),averaged.absolute_age(depth_index)],averaged.kim_oneil_temperature(depth_index)+[-2*averaged.kim_oneil_temperature_uncertainty(depth_index),2*averaged.kim_oneil_temperature_uncertainty(depth_index)],'Color',kim_oneil_colour);
    plot([averaged.absolute_age(depth_index),averaged.absolute_age(depth_index)],averaged.hansen_temperature(depth_index)+[-2*averaged.hansen_temperature_uncertainty(depth_index),2*averaged.hansen_temperature_uncertainty(depth_index)],'Color',hansen_colour);
    plot([averaged.absolute_age(depth_index),averaged.absolute_age(depth_index)],averaged.anderson_arthur_temperature(depth_index)+[-2*averaged.anderson_arthur_temperature_uncertainty(depth_index),2*averaged.anderson_arthur_temperature_uncertainty(depth_index)],'Color',anderson_arthur_colour);
end
oneil_handle = plot(averaged.absolute_age(10:end),averaged.oneil_temperature(10:end),'x-','Color',oneil_colour);
kim_oneil_handle = plot(averaged.absolute_age(10:end),averaged.kim_oneil_temperature(10:end),'x-','Color',kim_oneil_colour);
hansen_handle = plot(averaged.absolute_age(10:end),averaged.hansen_temperature(10:end),'x-','Color',hansen_colour);
anderson_arthur_handle = plot(averaged.absolute_age(10:end),averaged.anderson_arthur_temperature(10:end),'x-','Color',anderson_arthur_colour);


xlabel("Age");
ylabel("Temperature ^{\circ}C");

set(gca,"XDir","Reverse");

legend_handle = legend([oneil_handle,kim_oneil_handle,hansen_handle,anderson_arthur_handle],["O'Neil","Kim + O'Neil","Hansen","Anderson+Arthur"],"Location","NorthWest","Box","Off");

exportgraphics(gcf,"./../Plots/Age_Temperature_Calibrations.png");


figure(2);
clf
hold on

clf
hold on

plot(age_midpoints(1:9),oneil_difference(1:9),'.-','Color',oneil_colour);
plot(age_midpoints(1:9),kim_oneil_difference(1:9),'.-','Color',kim_oneil_colour);
plot(age_midpoints(1:9),hansen_difference(1:9),'.-','Color',hansen_colour);
plot(age_midpoints(1:9),anderson_arthur_difference(1:9),'.-','Color',anderson_arthur_colour);

plot(age_midpoints(9:end),oneil_difference(9:end),'x-','Color',oneil_colour);
plot(age_midpoints(9:end),kim_oneil_difference(9:end),'x-','Color',kim_oneil_colour);
plot(age_midpoints(9:end),hansen_difference(9:end),'x-','Color',hansen_colour);
plot(age_midpoints(9:end),anderson_arthur_difference(9:end),'.-','Color',anderson_arthur_colour);

xlabel("Age");
ylabel("\DeltaTemperature ^{\circ}C");

set(gca,"XDir","Reverse");

exportgraphics(gcf,"./../Plots/Age_deltaTemperature_Calibrations.png");