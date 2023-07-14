
clc;
clear all;

addpath(genpath(pwd)); % add current directory to path (mcmcstat is included)
addpath('./simulator/');
addpath('./utils/');
addpath('./mcmcstat/');


%% functions inclusion

%actual parameters used in the simulation
r = 0.2;
phi = 1.3e-8;
tau = 1.25;
beta = 200;
NE = 50;

%initial conditions
moi_mean = 0.1;
S0 = 1e8;
V0 = S0*moi_mean;
y(1,1) = S0;
y(1,2:NE+1) = 0;
y(1,NE+2) = 0 ;
y(1,NE+3) = V0;
time_free_phages = 0:0.1:2;



theta = [r,phi,tau,beta];

dilution_factor = 100;
[time,y_series_simulated,time_abs,pre_dil] = one_step_simulate_particular_points(time_free_phages,y,theta,NE, dilution_factor);

%%

figure(1)
plot(time,y_series_simulated,'o');
set(gca, 'YScale', 'log');
xlabel('time (hrs)');
ylabel('Free virions (/ml)')


%% Pick n points uniformly, N - n points remaining.
n = 10; % I choose 4 points to begin with and divide the time interval of experiment.
N = 15; % total points I will use for inference.    

t_selected = time(1:abs(length(time)/n):end);
V_selected = y_series_simulated_particular_points(NE+3,1:abs(length(time)/n):end);

% I will add noise later.

%% inference step 1 -- using n points.

data.ydata = V_selected;
data.xdata = t_selected;

model.ssfun = @(theta,data)  error_function_NE_varies(theta,data,S0,V0);

r_guess = 0.2;
phi_guess = 1.3e-8;
tau_guess = 1.2;
beta_guess = 200;
NE_guess = 50;

params = {
% initial values for the model states
    {'r', r_guess, 0, 0.5, r_guess,0.05}
    {'phi', phi_guess,  1e-10, 1e-6, phi_guess, 1e-7 }
    {'tau', tau_guess,   0.25, 5, tau_guess, 1 }
    {'beta', beta_guess, 0, 700, beta_guess, 100}
    {'NE',NE_guess,5,200,NE_guess,50};
    };

options.nsimu = 10000;
[results, chain, s2chain] = mcmcrun(model,data,params,options);

burn = options.nsimu*0.5;
theta_inferred = median(chain(burn:end,:));
NE_inferred = round(theta_inferred(5));

%% simulate again from inferred.

y0(1) = S0;
y0(2:NE_inferred+2) = 0;
y0(NE_inferred+3) = V0;

y00(1) = S0;
y00(2:52) = 0;
y00(53) = V0;



[time_inferred,y_series_inferred] = one_step_simulate(time_free_phages,y0,theta_inferred(1:4),NE_inferred,dilution_factor);
%[time_inferred,y_series_inferred] = one_step_simulate(time_free_phages,y00,[0.16 1.3e-8 0.25 102],17,dilution_factor);

%% plots

figure(2);
plot(time_inferred,y_series_inferred(NE_inferred+3,:),'-k');% hold on;plot(t_selected,V_selected,'r*')
set(gca, 'YScale', 'log');
xlabel('time (hrs)');
ylabel('Free virions (/ml)')
