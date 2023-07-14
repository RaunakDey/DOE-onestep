function ss = error_function_NE_varies(theta,data,S0,V0)
%Error function for mcmc

time_free_phages = data.xdata;
free_phages_mean = data.ydata;

dilution_factor = 100;

NE = round(theta(5));
%trial
[~,y_series2] = one_step_rescaled(time_free_phages,theta,NE,dilution_factor,S0,V0);
virus_simu = y_series2(end,:);

ss = sum((log(free_phages_mean) - log(y_series2(end,:))').^2);




end

