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
evolutions_preperturbation.subsample_boolean = logical(repmat([ones(9,1);zeros(size(evolutions.pH,1)-9,1)],1,size(evolutions.pH,2)));
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

co2_change_distribution = Geochemistry_Helpers.Distribution().fromSamples(-10:0.01:10,co2_doublings);
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
figure(1);
clf
hold on

plot([0,5],[0,2.5],'--','Color',[0.5,0.5,0.5]);
text(0.6,6.2,num2str(2.5/5)+" ^{\circ}C/doubling",'Rotation',75);
plot([0,5],[0,5],'--','Color',[0.5,0.5,0.5]);
text(1.4,6.3,num2str(5/5)+" ^{\circ}C/doubling",'Rotation',65);
plot([0,5],[0,10],'--','Color',[0.5,0.5,0.5]);
text(3.2,6.8,num2str(10/5)+" ^{\circ}C/doubling",'Rotation',50);
plot([0,5],[0,20],'--','Color',[0.5,0.5,0.5]);
text(4.2,4.5,num2str(20/5)+" ^{\circ}C/doubling",'Rotation',30);
plot([0,5],[0,40],'--','Color',[0.5,0.5,0.5]);
text(4.2,2.4,num2str(40/5)+" ^{\circ}C/doubling",'Rotation',20);

patch([co2_change_uncertainties(1),co2_change_uncertainties(1),co2_change_uncertainties(2),co2_change_uncertainties(2)],[temperature_change_uncertainties(1),temperature_change_uncertainties(2),temperature_change_uncertainties(2),temperature_change_uncertainties(1)],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.5);
text(2.25,4.2,"TJ Range",'Color',Geochemistry_Helpers.Colour.Colour("Black").darken(1.1).saturate(1.1).rgb,'FontWeight','Bold')

xlabel("\DeltaCO_2 (doublings)");
ylabel("\DeltaTemperature (^{\circ}C)");

xlim([0,5]);
ylim([0,8]);

exportgraphics(gcf,"./../../Figures/Climate_Sensitivity.png");
