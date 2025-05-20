clear

%% Load data
data_directory = "./../../../Data/";
boron_data = readtable(data_directory+"/Boron/TJ_d11B_d18O_d13C.xlsx","Sheet","Temperature_Calibrations");
good_boron_data = boron_data(~boron_data.diagenetic_alteration & ~boron_data.al_ca_reject,:);
rejected_boron_data = boron_data(boron_data.diagenetic_alteration | boron_data.al_ca_reject,:);

prior_file = jsondecode(fileread(data_directory+"/pH_change/prior.json"));
posterior_file = jsondecode(fileread(data_directory+"/pH_change/posterior.json"));
interpolation_ages = posterior_file.age;
prior = prior_file.prior;
posterior = posterior_file.posterior;

quantile_level = [2.5,97.5];
quantiles = (quantile_level/100);

d11B_evolutions_prior = [prior.d11B_measured];
for distribution_index = 1:size(d11B_evolutions_prior,1)
    evolutions.d11B_distributions_prior(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],d11B_evolutions_prior(distribution_index,:));
end

d11B_evolutions_posterior = [posterior.d11B_measured];

for distribution_index = 1:size(d11B_evolutions_posterior,1)
    evolutions.d11B_distributions_posterior(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],d11B_evolutions_posterior(distribution_index,:));
end
evolutions.d11B_mean_prior = evolutions.d11B_distributions_prior.mean();
evolutions.d11B_mean_posterior = evolutions.d11B_distributions_posterior.mean();

evolutions.d11B_quantiles_prior = evolutions.d11B_distributions_prior.quantile(quantiles);
evolutions.d11B_quantiles_posterior = evolutions.d11B_distributions_posterior.quantile(quantiles);

%%
clf
hold on
errorbar(good_boron_data.age,good_boron_data.d11B,good_boron_data.d11B_uncertainty,"LineStyle",'none',"Color",[0,0,0]);
scatter(good_boron_data.age,good_boron_data.d11B,15,"MarkerFaceColor",[0,0,0],"MarkerEdgeColor",[0,0,0])

errorbar(rejected_boron_data.age,rejected_boron_data.d11B,rejected_boron_data.d11B_uncertainty,"LineStyle",'none',"Color",[0,0,0]);
scatter(rejected_boron_data.age,rejected_boron_data.d11B,15,"MarkerFaceColor",[0,0,0],"MarkerEdgeColor",[0,0,0],"Marker","x")

p1 = patch([interpolation_ages;flipud(interpolation_ages)],[evolutions.d11B_quantiles_prior(:,1);flipud(evolutions.d11B_quantiles_prior(:,2))],"red","FaceAlpha",0.1);
plot(interpolation_ages,evolutions.d11B_mean_prior,"Color",[0.8,0.0,0.0],"LineWidth",4,"LineStyle",":");

p2 = patch([interpolation_ages;flipud(interpolation_ages)],[evolutions.d11B_quantiles_posterior(:,1);flipud(evolutions.d11B_quantiles_posterior(:,2))],"red","FaceAlpha",0.35);
plot(interpolation_ages,evolutions.d11B_mean_posterior,"Color",[0.7,0.0,0.0],"LineWidth",4);

% plot(interpolation_ages,d11B_evolutions_prior(:,1:5),"Color",[0.0,0.0,0.7]);
% plot(interpolation_ages,d11B_evolutions_posterior(:,1:5),"Color",[0.7,0.0,0.0]);

xlim([201.275,201.525]);
xlabel("Age (Ma)");
ylabel("\delta^{11}B (‰)");

current_axis = gca;
current_axis.YAxis.Exponent = 0;

xlim([min(interpolation_ages),max(interpolation_ages)]);
set(gca,'XDir','Reverse');

legend([p1,p2],["Prior","Posterior"]);