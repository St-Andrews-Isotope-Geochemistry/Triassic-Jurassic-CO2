clear

%% Load data
data_directory = "./../../../Data/";
boron_data = readtable(data_directory+"/Boron/TJ_d11B.xlsx","Sheet","Delta_Temperature");

prior_file = jsondecode(fileread(data_directory+"/pH_change/prior.json"));
posterior_file = jsondecode(fileread(data_directory+"/pH_change/posterior.json"));
interpolation_ages = posterior_file.age;
prior = prior_file.prior;
posterior = posterior_file.posterior;

quantile_level = [2.5,97.5];
quantiles = (quantile_level/100);

alkalinity_evolutions_prior = [prior.alkalinity];
for distribution_index = 1:size(alkalinity_evolutions_prior,1)
    evolutions.alkalinity_distributions_prior(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],alkalinity_evolutions_prior(distribution_index,:));
end

alkalinity_evolutions_posterior = [posterior.alkalinity];

for distribution_index = 1:size(alkalinity_evolutions_posterior,1)
    evolutions.alkalinity_distributions_posterior(distribution_index) = Geochemistry_Helpers.Distribution.fromSamples([],alkalinity_evolutions_posterior(distribution_index,:));
end
evolutions.alkalinity_mean_prior = evolutions.alkalinity_distributions_prior.mean();
evolutions.alkalinity_mean_posterior = evolutions.alkalinity_distributions_posterior.mean();

evolutions.alkalinity_quantiles_prior = evolutions.alkalinity_distributions_prior.quantile(quantiles);
evolutions.alkalinity_quantiles_posterior = evolutions.alkalinity_distributions_posterior.quantile(quantiles);

%%
clf
hold on
p1 = patch([interpolation_ages;flipud(interpolation_ages)],[evolutions.alkalinity_quantiles_prior(:,1);flipud(evolutions.alkalinity_quantiles_prior(:,2))],"red","FaceAlpha",0.1);
plot(interpolation_ages,evolutions.alkalinity_mean_prior,"Color",[0.8,0.0,0.0],"LineWidth",4,"LineStyle",":");

p2 = patch([interpolation_ages;flipud(interpolation_ages)],[evolutions.alkalinity_quantiles_posterior(:,1);flipud(evolutions.alkalinity_quantiles_posterior(:,2))],"red","FaceAlpha",0.35);
plot(interpolation_ages,evolutions.alkalinity_mean_posterior,"Color",[0.7,0.0,0.0],"LineWidth",4);

% plot(interpolation_ages,alkalinity_evolutions_prior(:,1:5),"Color",[0.0,0.0,0.7]);
% plot(interpolation_ages,alkalinity_evolutions_posterior(:,1:5),"Color",[0.7,0.0,0.0]);

xlim([201.275,201.525]);
xlabel("Age (Ma)");
ylabel("Alkalinity (\mumol/kg)");

current_axis = gca;
current_axis.YAxis.Exponent = 0;

xlim([min(interpolation_ages),max(interpolation_ages)]);
set(gca,'XDir','Reverse');

legend([p1,p2],["Prior","Posterior"]);