clear

%% Load data
d18O_d13C = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Delta_Temperature");
d18O_d13C_averaged = readtable("./../../Data/TJ_d18O_d13C.xlsx","Sheet","Averaged");

boron_data = readtable("./../../Data/TJ_d11B.xlsx","Sheet","Delta_Temperature");
evolutions = getShapedEvolutions("./../../Data/TJ_CO2_Evolutions.csv");

interpolation_ages = unique(sort([boron_data.age',linspace(min(boron_data.age),max(boron_data.age),80)]));
evolutions.age = repmat(interpolation_ages',1,size(evolutions.pH,2));

%% Get the initial subsample
evolutions_preperturbation.subsample_boolean = evolutions.age>=boron_data.age(9);
evolutions_preperturbation.pH = reshape(evolutions.pH(evolutions_preperturbation.subsample_boolean),9,[]);
evolutions_preperturbation.co2 =  reshape(evolutions.co2(evolutions_preperturbation.subsample_boolean),9,[]);

evolutions_preperturbation.co2_mean = mean(evolutions_preperturbation.co2);
evolutions_preperturbation.co2_initial = evolutions_preperturbation.co2(1,:);

d11B_minimum = logical([zeros(9,1);boron_data.d11B(10:end)==min(boron_data.d11B(10:end))]);
evolutions_d11B_minimum.subsample_boolean = repmat(d11B_minimum,1,size(evolutions.pH,2));
evolutions_d11B_minimum.pH = reshape(evolutions.pH(evolutions_d11B_minimum.subsample_boolean),1,[]);
evolutions_d11B_minimum.co2 =  reshape(evolutions.co2(evolutions_d11B_minimum.subsample_boolean),1,[]);

co2_change = evolutions_d11B_minimum.co2-evolutions_preperturbation.co2_initial;
co2_fractional_change = evolutions_d11B_minimum.co2./evolutions_preperturbation.co2_initial;
co2_doublings = log2(co2_fractional_change);
co2_doublings = co2_doublings(imag(co2_doublings)==0);

co2_change_distribution = Geochemistry_Helpers.Distribution().fromSamples(-10:0.01:10,co2_doublings).normalise();
co2_change_sampler = Geochemistry_Helpers.Sampler(co2_change_distribution,"latin_hypercube");

co2_change_sampler.getSamples(numel(co2_change)).shuffle();

%% Generate random possible temperature changes
temperature_change_parameters = [4.2,1.1];
temperature_change_distribution = Geochemistry_Helpers.Distribution(0:0.1:8,"Gaussian",temperature_change_parameters).normalise();
temperature_change_sampler = Geochemistry_Helpers.Sampler(temperature_change_distribution,"latin_hypercube");

temperature_change_sampler.getSamples(numel(co2_change)).shuffle();

%%
co2_change_median = co2_change_distribution.median();
co2_change_uncertainties = [co2_change_distribution.quantile(0.025),co2_change_distribution.quantile(0.975)];

temperature_change_median = temperature_change_distribution.median();
temperature_change_uncertainties = [temperature_change_distribution.quantile(0.025),temperature_change_distribution.quantile(0.975)];

%%
climate_sensivity = temperature_change_sampler.samples./co2_change_sampler.samples;
climate_sensitivity_distribution = Geochemistry_Helpers.Distribution().fromSamples(0:0.1:20,climate_sensivity).normalise();

temperature_change_axis_matrix = repmat(temperature_change_distribution.bin_midpoints,numel(co2_change_distribution.probabilities),1);
temperature_change_probability_matrix = repmat(temperature_change_distribution.probabilities,numel(co2_change_distribution.probabilities),1);

co2_change_axis_matrix = repmat(co2_change_distribution.bin_midpoints,1,numel(temperature_change_distribution.probabilities));
co2_change_probability_matrix = repmat(co2_change_distribution.probabilities,1,numel(temperature_change_distribution.probabilities));

climate_sensitivity_joint = temperature_change_probability_matrix.*co2_change_probability_matrix;
filter = ones(20,20);
normalised_filter = filter./sum(filter);
climate_sensitivity_joint_smooth = conv2(climate_sensitivity_joint,normalised_filter,'same');
climate_sensitivity_joint_smooth_normalised = climate_sensitivity_joint_smooth./sum(climate_sensitivity_joint_smooth(:));

%%
colourmap = Geochemistry_Helpers.Colour.Map("b",[Geochemistry_Helpers.Colour.Colour("Crimson","rgb",2),...
                                                 Geochemistry_Helpers.Colour.Colour("Crimson","rgb",1.8),...
                                                 Geochemistry_Helpers.Colour.Colour("Orange","rgb",1.2),...
                                                 Geochemistry_Helpers.Colour.Colour("Khaki","rgb",0.5),...
                                                 Geochemistry_Helpers.Colour.Colour("White","rgb",0)]);
expanded_colourmap = colourmap.getColours(100);

%% 
figure(1);
clf
hold on

pcolor_handle = pcolor(co2_change_axis_matrix,temperature_change_axis_matrix,climate_sensitivity_joint_smooth_normalised);
colormap(expanded_colourmap.colours.rgb);
set(pcolor_handle,'EdgeColor','None');
shading interp;

xlabel("\DeltaCO_2 (doublings)");
ylabel("\DeltaTemperature (^{\circ}C)");

plot([0,5],[0,2.5],'--','Color',[0.5,0.5,0.5]);
text(4.2,2.4,num2str(2.5/5)+" ^{\circ}C/doubling",'Rotation',20);
plot([0,5],[0,5],'--','Color',[0.5,0.5,0.5]);
text(4.2,4.5,num2str(5/5)+" ^{\circ}C/doubling",'Rotation',30);
plot([0,5],[0,10],'--','Color',[0.5,0.5,0.5]);
text(3.2,6.8,num2str(10/5)+" ^{\circ}C/doubling",'Rotation',50);
plot([0,5],[0,20],'--','Color',[0.5,0.5,0.5]);
text(1.4,6.3,num2str(20/5)+" ^{\circ}C/doubling",'Rotation',65);
plot([0,5],[0,40],'--','Color',[0.5,0.5,0.5]);
text(0.6,6.2,num2str(40/5)+" ^{\circ}C/doubling",'Rotation',75);


xlim([0,5]);
ylim([0,8]);

exportgraphics(gcf,"./../../Figures/Climate_Sensitivity_Probabilistic.png");
