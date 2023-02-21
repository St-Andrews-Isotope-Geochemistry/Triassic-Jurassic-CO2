clear

starting_pH = 7.5:0.01:9;
starting_co2s = BuCC.d11BCO2().create(numel(starting_pH));
starting_co2s.collate("boron").assignToAll("d11B_4",Geochemistry_Helpers.delta("B",NaN));
starting_co2s.collate("boron").assignToAll("d11B_sw",Geochemistry_Helpers.delta("B",33.32));
starting_co2s.collate("boron").collate("pH").assignToEach("pValue",starting_pH);
starting_co2s.collate("boron").assignToAll("epsilon",27.2);

starting_co2s.collate("carbonate_chemistry").assignToAll("temperature",21.4);
starting_co2s.collate("carbonate_chemistry").assignToAll("salinity",35);
starting_co2s.collate("carbonate_chemistry").assignToAll("oceanic_pressure",0);
starting_co2s.collate("carbonate_chemistry").assignToAll("atmospheric_pressure",1);
starting_co2s.collate("carbonate_chemistry").assignToAll("calcium",17);
starting_co2s.collate("carbonate_chemistry").assignToAll("magnesium",28);

starting_co2s.collate("carbonate_chemistry").collate("atmospheric_co2").assignToAll("partial_pressure",500e-6);

starting_co2s.calculate();
starting_co2 = starting_co2s.collate("carbonate_chemistry").collate("atmospheric_co2").collate("mole_fraction");

starting_d11B_4 = starting_co2s.collate("boron").collate("d11B_4").collate("value");

d11B_4_change = [2,4,8,3.5126];
final_d11B4s = Geochemistry_Helpers.delta().create([numel(starting_pH),numel(d11B_4_change)]);
final_d11B4s.assignToAll("standard","B");
final_d11B4s.assignToEach("value",starting_d11B_4-d11B_4_change);

final_co2s = BuCC.d11BCO2().create([numel(starting_pH),numel(d11B_4_change)]);
final_co2s.collate("boron").assignToEach("d11B_4",final_d11B4s);
final_co2s.collate("boron").collate("d11B_sw").assignToAll("value",33.32);
final_co2s.collate("boron").assignToAll("epsilon",27.2);
final_co2s.collate("boron").assignToEach("pKb",repmat(starting_co2s.collate("boron").collate("pKb"),1,numel(d11B_4_change)));
final_co2s.collate("boron").collate("pH").assignToAll("pValue",NaN);

final_co2s.collate("carbonate_chemistry").assignToAll("temperature",21.4);
final_co2s.collate("carbonate_chemistry").assignToAll("salinity",35);
final_co2s.collate("carbonate_chemistry").assignToAll("oceanic_pressure",0);
final_co2s.collate("carbonate_chemistry").assignToAll("atmospheric_pressure",1);
final_co2s.collate("carbonate_chemistry").assignToAll("calcium",17);
final_co2s.collate("carbonate_chemistry").assignToAll("magnesium",28);

final_co2s.collate("carbonate_chemistry").assignToEach("alkalinity",repmat(starting_co2s.collate("carbonate_chemistry").collate("alkalinity")',1,numel(d11B_4_change)));
final_co2s.collate("carbonate_chemistry").assignToAll("units"," mol/kg");

final_co2s.calculate();

%
final_pH = final_co2s.collate("carbonate_chemistry").collate("pH").collate("pValue");
real_pH = imag(final_pH)==0;
final_pH(~real_pH) = NaN;

final_co2 = final_co2s.collate("carbonate_chemistry").collate("atmospheric_co2").collate("mole_fraction");
final_co2(~real_pH) = NaN;

pH_change = (final_pH-starting_pH');
co2_change = (final_co2-starting_co2)*1e6;

%%
pH_colours{3} = [0.8,0,0];
pH_colours{2} = [0.2,0.6,0.8];
pH_colours{1} = [0.3,0.7,0.3];

clf
subplot_handles(1) = subplot(2,1,1);
hold on
for change_index = 1:size(pH_change,2)-1
    plot(starting_pH,pH_change(:,change_index),'Color',pH_colours{change_index},'LineWidth',2);
end
plot(starting_pH,pH_change(:,end),'k:','LineWidth',2);

set(gca,'XTick',[]);
ylabel("\DeltapH");

xlim([min(starting_pH),max(starting_pH)]);
ylim([-1.5,0]);

legend_handle = legend(["2 ‰","4 ‰","8 ‰","3.51 ‰"],'Box','Off','Location','SouthEast');
title(legend_handle,'\delta^{11}B_4 change');

subplot_handles(2) = subplot(2,1,2);
hold on
for change_index = 1:size(pH_change,2)-1
    plot(starting_pH,co2_change(:,change_index),'Color',pH_colours{change_index},'LineWidth',2);
end
plot(starting_pH,co2_change(:,end),'k:','LineWidth',2);

xlabel("Initial pH");
ylabel("\DeltaxCO_2 (ppm)");

co2_axis = gca;
co2_axis.YAxis.Exponent = 0;

ylim([0,10000]);

text(8.75,9000,char(949)+" = 27.2"+char(8240));
text(8.63,7500,"\delta^{11}B_{sw} = 33.2"+char(8240));

linkaxes(subplot_handles,"x");

% exportgraphics(gcf,"./../../Figures/Initial_pH_Delta_pH_Delta_CO2.png");

