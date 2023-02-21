clear
data_directory = "./../../../Data/";

minimum_pH_ensemble = readmatrix(data_directory+"/Minimum_pH_Change/TJ_Minimum_pH_Ensemble.csv");
pH = minimum_pH_ensemble(3:4,:);

pH_initial = pH(1,:);
pH_after = pH(2,:);

% If a value in either array is imaginary then we don't want it
real_pH = imag(pH_after)==0 & imag(pH_initial)==0;

% Create distributions from the samples
pH_initial_distribution = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,pH_initial(real_pH)).normalise();
pH_after_distribution = Geochemistry_Helpers.Distribution.fromSamples(6:0.01:9,pH_after(real_pH)).normalise();

pH_change = pH_after(real_pH) - pH_initial(real_pH);
pH_change_distribution = Geochemistry_Helpers.Distribution.fromSamples(-1:0.01:1,pH_change).normalise();

%% Statistics
% Can also calculate metrics to describe the distribution
pH_change_99 = pH_change_distribution.quantile(0.99);
pH_change_975 = pH_change_distribution.quantile(0.975);
pH_change_95 = pH_change_distribution.quantile(0.95);
pH_change_50 = pH_change_distribution.quantile(0.5);
pH_change_025 = pH_change_distribution.quantile(0.025);

pH_initial_99 = pH_initial_distribution.quantile(0.99);
pH_initial_975 = pH_initial_distribution.quantile(0.975);
pH_initial_95 = pH_initial_distribution.quantile(0.95);
pH_initial_50 = pH_initial_distribution.quantile(0.5);
pH_initial_025 = pH_initial_distribution.quantile(0.025);

pH_after_99 = pH_after_distribution.quantile(0.99);
pH_after_975 = pH_after_distribution.quantile(0.975);
pH_after_95 = pH_after_distribution.quantile(0.95);
pH_after_50 = pH_after_distribution.quantile(0.5);
pH_after_025 = pH_after_distribution.quantile(0.025);

%%
figure(1)
clf
pH_change_distribution.plot();

xlim([-0.5,-0.25]);
xlabel("\DeltapH");
ylabel("Probability");


%%
filename = "/Minimum_pH_Change/Metrics.json";

fileID = fopen(data_directory+filename,"w");
fwrite(fileID,"{"+newline+string(char(9))+'"initial":');
fwrite(fileID,"{"+newline+string(char(9))+string(char(9))+'"pH_025":'+pH_initial_025+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_median":'+pH_initial_50+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_975":'+pH_initial_975+newline+string(char(9))+"},"+newline);
fwrite(fileID,string(char(9))+'"after":');
fwrite(fileID,"{"+newline+string(char(9))+string(char(9))+'"pH_025":'+pH_after_025+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_median":'+pH_after_50+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_975":'+pH_after_975+newline+string(char(9))+"},"+newline);
fwrite(fileID,string(char(9))+'"change":');
fwrite(fileID,"{"+newline+string(char(9))+string(char(9))+'"pH_025":'+pH_change_025+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_median":'+pH_change_50+",");
fwrite(fileID,newline+string(char(9))+string(char(9))+'"pH_95":'+pH_change_975+newline+string(char(9))+"}"+newline);
fwrite(fileID,"}");
fclose(fileID);

