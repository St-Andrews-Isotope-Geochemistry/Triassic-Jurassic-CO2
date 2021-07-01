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
evolutions.index = repmat((1:size(evolutions.pH,1))',1,size(evolutions.pH,2));

clear raw_evolutions reshaped_evolutions

%% Create colours
colours = Geochemistry_Helpers.Colour.Map("rb",[Geochemistry_Helpers.Colour.Colour("Firebrick","ryb",0),...
                                                Geochemistry_Helpers.Colour.Colour("Orchid","ryb",0.5),...
                                                Geochemistry_Helpers.Colour.Colour("Steelblue","ryb",1)]);

%%
evolutions_preperturbation.subsample_boolean = evolutions.index<=9;
evolutions_preperturbation.pH = reshape(evolutions.pH(evolutions_preperturbation.subsample_boolean),9,[]);
evolutions_preperturbation.co2 =  reshape(evolutions.co2(evolutions_preperturbation.subsample_boolean),9,[]);

evolutions_preperturbation.pH_mean = mean(evolutions_preperturbation.pH);
evolutions_preperturbation.co2_mean = mean(evolutions_preperturbation.co2);

evolutions_difference.pH = evolutions.pH-repmat(evolutions_preperturbation.pH_mean,22,1);
evolutions_difference.co2 = evolutions.co2-repmat(evolutions_preperturbation.co2_mean,22,1);

%%
pH_bin_width = 0.5;
pH_bin_edges = 7.6:pH_bin_width:8.6-pH_bin_width;

colourmap = colours.getColours(numel(pH_bin_edges));
% colourmap_2 = colours.colours(1).makePalette("triad").getColours(numel(pH_bin_edges));

for pH_bin_index = 1:numel(pH_bin_edges)
    evolutions_split(pH_bin_index).subsample_boolean = repmat(evolutions_preperturbation.pH_mean>=pH_bin_edges(pH_bin_index) & evolutions_preperturbation.pH_mean< pH_bin_edges(pH_bin_index)+pH_bin_width,size(evolutions.pH,1),1);
    evolutions_split(pH_bin_index).pH_difference = reshape(evolutions_difference.pH(evolutions_split(pH_bin_index).subsample_boolean),22,[]);
    evolutions_split(pH_bin_index).co2_difference =  reshape(evolutions_difference.co2(evolutions_split(pH_bin_index).subsample_boolean),22,[]);
    
    for age_index = 1:numel(boron_data.age)-1
        evolutions_split(pH_bin_index).delta_pH_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(-1:0.01:1,evolutions_split(pH_bin_index).pH_difference(age_index,:)).normalise();
        evolutions_split(pH_bin_index).delta_co2_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(-10000:100:10000,evolutions_split(pH_bin_index).co2_difference(age_index,:)).normalise();
    end
    evolutions_split(pH_bin_index).pH_difference_median = evolutions_split(pH_bin_index).delta_pH_distributions.median();
    evolutions_split(pH_bin_index).pH_difference_95 = [evolutions_split(pH_bin_index).delta_pH_distributions.quantile(0.025),evolutions_split(pH_bin_index).delta_pH_distributions.quantile(0.975)];
end


%% 
figure(1);
clf
hold on
for evolution_index = 1:numel(evolutions_split)
    patch([boron_data.age(1:end-1);flipud(boron_data.age(1:end-1))],[evolutions_split(evolution_index).pH_difference_95(:,1);flipud(evolutions_split(evolution_index).pH_difference_95(:,2))],colourmap.colours(evolution_index).rgb,'EdgeColor','None','FaceAlpha',0.5);
    plot_handles(evolution_index) = plot(boron_data.age(1:end-1),evolutions_split(evolution_index).pH_difference_median,'Color',colourmap.colours(evolution_index).rgb);
end

set(gca,'XDir','Reverse','XTick',[]);
ylabel("\DeltapH");

legend_handle = legend(plot_handles,string(num2str(pH_bin_edges'+pH_bin_width/2)),'Location','SouthWest');
title(legend_handle,"Initial pH");

exportgraphics(gcf,"./../../Figures/Delta_pH_Split_Initial_pH.png");

%%
pH_bin_width = 0.2;
pH_bin_edges = 7.6:pH_bin_width:8.4-pH_bin_width;

colourmap = colours.getColours(numel(pH_bin_edges));
% colourmap_2 = colours.colours(1).makePalette("triad").getColours(numel(pH_bin_edges));

for pH_bin_index = 1:numel(pH_bin_edges)
    evolutions_split(pH_bin_index).subsample_boolean = repmat(evolutions_preperturbation.pH_mean>=pH_bin_edges(pH_bin_index) & evolutions_preperturbation.pH_mean< pH_bin_edges(pH_bin_index)+pH_bin_width,size(evolutions.pH,1),1);
    evolutions_split(pH_bin_index).pH_difference = reshape(evolutions_difference.pH(evolutions_split(pH_bin_index).subsample_boolean),22,[]);
    evolutions_split(pH_bin_index).co2_difference =  reshape(evolutions_difference.co2(evolutions_split(pH_bin_index).subsample_boolean),22,[]);
    
    for age_index = 1:numel(boron_data.age)-1
        evolutions_split(pH_bin_index).delta_pH_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(-1:0.01:1,evolutions_split(pH_bin_index).pH_difference(age_index,:)).normalise();
        evolutions_split(pH_bin_index).delta_co2_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(-10000:100:10000,evolutions_split(pH_bin_index).co2_difference(age_index,:)).normalise();
    end
    evolutions_split(pH_bin_index).pH_difference_median = evolutions_split(pH_bin_index).delta_pH_distributions.median();
    evolutions_split(pH_bin_index).pH_difference_95 = [evolutions_split(pH_bin_index).delta_pH_distributions.quantile(0.025),evolutions_split(pH_bin_index).delta_pH_distributions.quantile(0.975)];
end

figure(2);
clf
hold on
for evolution_index = 1:numel(evolutions_split)
    plot_handles(evolution_index) = plot(boron_data.age(1:end-1),evolutions_split(evolution_index).pH_difference_median,'Color',colourmap.colours(evolution_index).rgb,'LineWidth',2);
end

set(gca,'XDir','Reverse','XTick',[]);
ylabel("\DeltapH");

legend_handle = legend(plot_handles,string(num2str(pH_bin_edges'+pH_bin_width/2)),'Location','SouthWest');
title(legend_handle,"Initial pH");

exportgraphics(gcf,"./../../Figures/Delta_pH_Split_Initial_pH_Window.png");