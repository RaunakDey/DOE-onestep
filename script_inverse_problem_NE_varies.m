clc;
clear all;
addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('/Users/rdey33/Downloads/MATLAB_DRIVE/mcmcstat/mcmcstat');

name = "CBA18-3_18";
load('./data/'+name+'.mat');
[indi_i,indi_j]= name_matrix_function(name);
time_free_phages = time_free_phages/60;

load('parameters.mat','pars');


%% collecting statistics

free_phages_mean = mean(free_phages,2);
seed = 7;
filestr = "one_step_results/"+name+"/local/seed"+string(seed);

data.ydata = free_phages_mean;
data.xdata = time_free_phages;

%% functions inclusion
%NE =85; %number of exposed class

moi_mean = mean(moi);
S0 = 1e8;
V0 = S0*moi_mean;



% theta given by OSU lab
beta_osu = prior_values.burst_size;
r_osu = pars.r(indi_i);
phi_osu = pars.phi(indi_i,indi_j);
tau_osu = prior_values.latent_period/60  ;

theta_osu = [r_osu,phi_osu,tau_osu,beta_osu];


%% inference mcmc

model.ssfun = @(theta,data)  error_function_NE_varies(theta,data,S0,V0);


params = {
% initial values for the model states
    {'r', r_osu, 0, 0.5, r_osu,0.05}
    {'phi', phi_osu,  1e-10, 1e-6, phi_osu, 1e-7 }
    {'tau', tau_osu,   0.25, 5, tau_osu, 1 }
    {'beta', beta_osu, 0, 700, beta_osu, 100}
    {'NE',70,5,200,70,50};
    };



options.nsimu = 10000;
[results, chain, s2chain] = mcmcrun(model,data,params,options);

chainstats(chain,results);


%% plots
burn = options.nsimu*0.8;
format short E;
theta_inferred = median(chain(burn:end,:));

NE = round(theta_inferred(5));
y0(1) = S0;
y0(2:NE+2) = 0;
y0(NE+3) = V0;



%theta_osu(1) =  0.21;

dilution_factor = 100;
[time,y_series_inferred] = one_step_simulate(time_free_phages,y0,theta_inferred(1:4),NE,dilution_factor);
[time2,y_series_osu] = one_step_simulate(time_free_phages,y0,theta_osu,NE,dilution_factor);
%[~,virus_all] = confidence_interval(time_free_phages,data,chain(burn:end,:),dilution_factor);

%% figures

figure(2)
[~,something] = size(free_phages);
if something == 1
plot(time_free_phages,free_phages,'o','LineWidth', 2);
else
errorbar(time_free_phages,mean(free_phages,2),std(free_phages'),'LineWidth',2,'LineStyle','none','Marker','o');
end
xlabel('time (hours)','interpreter','none')
ylabel("Free phages (" +string(cell2mat(labels.units(1,3)))+")",'interpreter','none');
title(string(labels.phage)+'-'+string(labels.host));
hold on;
plot(time,y_series_inferred(end,:),'-k','LineWidth',2);hold on;
plot(time2,y_series_osu(end,:),'-r','LineWidth',2);
legend('One step data','SEIV with inferred parameters','SEIV with OSU parameters','Location','northoutside');
set(gca, 'YScale', 'log');
set(gca,'FontSize',18)
str = {['OSU'],['r = ',num2str(theta_osu(1)),' cells/hr'], ['\phi = ',num2str(theta_osu(2))], ['\tau = ',num2str(theta_osu(3)),' hr'], ['\beta = ',num2str(theta_osu(4))]};
annotation('textbox', [0.7, 0.45, 0.1, 0.1], 'String', str,'FontSize',11,'FitBoxToText','on');
str2 = {['Inferred'],['r = ',num2str(theta_inferred(1)),' cells/hr'], ['\phi = ',num2str(theta_inferred(2))], ['\tau = ',num2str(theta_inferred(3)),' hr'], ['\beta = ',num2str(theta_inferred(4))] ['N_E = ',num2str(round(theta_inferred(5)))]};
annotation('textbox', [0.7, 0.25, 0.1, 0.1], 'String', str2,'FontSize',11,'FitBoxToText','on');

saveas(gcf,'./results_same_phi/'+name+'_seed'+string(seed)+'_timeseries.png')





figure(3)
subplot(5,2,1)
plot(chain(:,1));xlabel('MCMC step');ylabel('r (cells/hr)');
subplot(5,2,2)
histogram(chain(burn:end,1),'Normalization','probability','DisplayStyle','stairs','NumBins',50); xlabel('r (cells/hr') ; ylabel('Counts');

subplot(5,2,3)
plot(log(chain(:,2))./log(10));xlabel('MCMC step');ylabel('\phi');
subplot(5,2,4)
histogram(log(chain(burn:end,2))./log(10),'Normalization','probability','DisplayStyle','stairs','NumBins',50); xlabel('\phi') ; ylabel('\phi');

subplot(5,2,5)
plot(chain(:,3));xlabel('MCMC step');ylabel('\tau (hr)');
subplot(5,2,6)
histogram(chain(burn:end,3),'Normalization','probability','DisplayStyle','stairs','NumBins',50); xlabel('\tau (hr)') ; ylabel('Counts');

subplot(5,2,7)
plot(chain(:,4));xlabel('MCMC step');ylabel('\beta');
subplot(5,2,8)
histogram(chain(burn:end,4),'Normalization','probability','DisplayStyle','stairs','NumBins',50); xlabel('\beta') ; ylabel('Counts');

subplot(5,2,9)
plot(chain(:,5));xlabel('MCMC step');ylabel('N_E');
subplot(5,2,10)
histogram(chain(burn:end,5),'Normalization','probability','DisplayStyle','stairs','NumBins',30); xlabel('N_E') ; ylabel('Counts');


saveas(gcf,'./results_same_phi/'+name+'_seed'+string(seed)+'_estimates.png')



save('./results_same_phi/'+name+'_'+seed+'-inferred.mat');

%% error function

for i=1:200
    error_stored(i) = error_function_NE_varies(chain(9000+i*5,:),data,S0,V0);
end

figure(4)
plot(9000+(1:200)*5,error_stored);xlabel('chain position');ylabel('error function');
