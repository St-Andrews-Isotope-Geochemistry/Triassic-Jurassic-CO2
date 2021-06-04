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
evolutions_subsample{1,1}.subsample_boolean = repmat(evolutions.co2(1,:)<=2000 & evolutions.saturation_state(1,:)<=8,size(evolutions.pH,1),1);
evolutions_subsample{1,1} = expandSubsampleBoolean(evolutions,evolutions_subsample{1,1});

evolutions_subsample{2,1}.subsample_boolean = repmat(evolutions.co2(1,:)>2000 & evolutions.saturation_state(1,:)<=8,size(evolutions.pH,1),1);
evolutions_subsample{2,1} = expandSubsampleBoolean(evolutions,evolutions_subsample{2,1});

evolutions_subsample{1,2}.subsample_boolean = repmat(evolutions.co2(1,:)<=2000 & evolutions.saturation_state(1,:)>8 & evolutions.saturation_state(1,:)<=10,size(evolutions.pH,1),1);
evolutions_subsample{1,2} = expandSubsampleBoolean(evolutions,evolutions_subsample{1,2});

evolutions_subsample{2,2}.subsample_boolean = repmat(evolutions.co2(1,:)>2000 & evolutions.saturation_state(1,:)>8 & evolutions.saturation_state(1,:)<=10,size(evolutions.pH,1),1);
evolutions_subsample{2,2} = expandSubsampleBoolean(evolutions,evolutions_subsample{2,2});

evolutions_subsample{1,3}.subsample_boolean = repmat(evolutions.co2(1,:)<=2000 & evolutions.saturation_state(1,:)>10,size(evolutions.pH,1),1);
evolutions_subsample{1,3} = expandSubsampleBoolean(evolutions,evolutions_subsample{1,3});

evolutions_subsample{2,3}.subsample_boolean = repmat(evolutions.co2(1,:)>2000 & evolutions.saturation_state(1,:)>10,size(evolutions.pH,1),1);
evolutions_subsample{2,3} = expandSubsampleBoolean(evolutions,evolutions_subsample{2,3});

%%
figure(1);
clf

subplot_handles(1) = subplot(2,3,1);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{1,1}.pH_uncertainties(:,1);flipud(evolutions_subsample{1,1}.pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{1,1}.pH_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(2) = subplot(2,3,2);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{1,2}.pH_uncertainties(:,1);flipud(evolutions_subsample{1,2}.pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{1,2}.pH_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(3) = subplot(2,3,3);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{1,3}.pH_uncertainties(:,1);flipud(evolutions_subsample{1,3}.pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{1,3}.pH_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(4) = subplot(2,3,4);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{2,1}.pH_uncertainties(:,1);flipud(evolutions_subsample{2,1}.pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{2,1}.pH_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(5) = subplot(2,3,5);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{2,2}.pH_uncertainties(:,1);flipud(evolutions_subsample{2,2}.pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{2,2}.pH_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(6) = subplot(2,3,6);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{2,3}.pH_uncertainties(:,1);flipud(evolutions_subsample{2,3}.pH_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{2,3}.pH_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);


linkaxes(subplot_handles,"y");
set(subplot_handles,'XDir','Reverse');
set(subplot_handles([2,3,5,6]),"ytick",[]);
set(subplot_handles([1,2,3]),"xtick",[]);

xlabel(subplot_handles([4,5,6]),"Age (Ma)");
ylabel(subplot_handles([1,4]),"pH");

text(subplot_handles(1),201.52,8.7,"Saturation state <8");
text(subplot_handles(2),201.57,8.7,"Saturation state >8 and <12");
text(subplot_handles(3),201.52,8.7,"Saturation state >12");

text(subplot_handles(1),201.63,7,"CO_2 <2000",'Rotation',90);
text(subplot_handles(4),201.63,7,"CO_2 >2000",'Rotation',90);

figure(2);
clf

subplot_handles(1) = subplot(2,3,1);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{1,1}.co2_uncertainties(:,1);flipud(evolutions_subsample{1,1}.co2_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{1,1}.co2_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(2) = subplot(2,3,2);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{1,2}.co2_uncertainties(:,1);flipud(evolutions_subsample{1,2}.co2_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{1,2}.co2_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(3) = subplot(2,3,3);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{1,3}.co2_uncertainties(:,1);flipud(evolutions_subsample{1,3}.co2_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{1,3}.co2_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(4) = subplot(2,3,4);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{2,1}.co2_uncertainties(:,1);flipud(evolutions_subsample{2,1}.co2_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{2,1}.co2_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(5) = subplot(2,3,5);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{2,2}.co2_uncertainties(:,1);flipud(evolutions_subsample{2,2}.co2_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{2,2}.co2_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);

subplot_handles(6) = subplot(2,3,6);
hold on
patch([boron_data.age;flipud(boron_data.age)],[evolutions_subsample{2,3}.co2_uncertainties(:,1);flipud(evolutions_subsample{2,3}.co2_uncertainties(:,2))],Geochemistry_Helpers.Colour.Colour("SteelBlue").rgb,'FaceAlpha',0.3,'EdgeColor',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);
plot(boron_data.age,evolutions_subsample{2,3}.co2_median,'Color',Geochemistry_Helpers.Colour.Colour("SteelBlue").darken(1.2).rgb);


linkaxes(subplot_handles,"y");
set(subplot_handles,'XDir','Reverse');
set(subplot_handles([2,3,5,6]),"ytick",[]);
set(subplot_handles([1,2,3]),"xtick",[]);

xlabel(subplot_handles([4,5,6]),"Age (Ma)");
ylabel(subplot_handles([1,4]),"CO_2");

text(subplot_handles(1),201.52,11000,"Saturation state <8");
text(subplot_handles(2),201.57,11000,"Saturation state >8 and <12");
text(subplot_handles(3),201.52,11000,"Saturation state >12");

text(subplot_handles(1),201.67,3000,"CO_2 <2000",'Rotation',90);
text(subplot_handles(4),201.67,3000,"CO_2 >2000",'Rotation',90);

%%
figure(1);
exportgraphics(gcf,"./../../Figures/Initial_CO2_SaturationState_pH.png");

figure(2);
exportgraphics(gcf,"./../../Figures/Initial_CO2_SaturationState_CO2.png");

%%
function subsample = expandSubsampleBoolean(evolutions,subsample)
    subsample.pH = reshape(evolutions.pH(subsample.subsample_boolean),22,[]);
    subsample.co2 =  reshape(evolutions.co2(subsample.subsample_boolean),22,[]);

    subsample.pH_distributions = Geochemistry_Helpers.Distribution().create(size(subsample.pH,1));
    subsample.co2_distributions = Geochemistry_Helpers.Distribution().create(size(subsample.pH,1));
    for age_index = 1:size(subsample.pH,1)
        subsample.pH_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(5:0.01:10,subsample.pH(age_index,:)).normalise();
        subsample.co2_distributions(age_index) = Geochemistry_Helpers.Distribution.fromSamples(0:100:10000,subsample.co2(age_index,:)).normalise();
    end

    subsample.pH_median = subsample.pH_distributions.median();
    subsample.pH_uncertainties = [subsample.pH_distributions.quantile(0.025),subsample.pH_distributions.quantile(0.975)];

    subsample.co2_median = subsample.co2_distributions.median();
    subsample.co2_uncertainties = [subsample.co2_distributions.quantile(0.025),subsample.co2_distributions.quantile(0.975)];
end