function [time,y] = one_step_rescaled(time_free_phages,theta,NE,dilution_factor,S0,V0)
%ONE_STEP_RESCALES Summary of this function goes here
%   Detailed explanation goes here


opts=odeset('AbsTol',1E-8,'RelTol',1E-8,'NonNegative',1,'MaxStep',0.25);

dt_dilution = 1/60;  
time_for_absorbtion = 0:dt_dilution:15/60;  %15 mins given to us

%V0 - data.ydata(1)

y0(1) = S0;
y0(2:NE+1) = 0;
y0(NE+2) = 0;
y0(NE+3) = V0;

solved = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), time_for_absorbtion',y0,opts);

%raunak's idea
after_dilution = solved.y(:,end)/dilution_factor;
%after_dilution(2:NE+2) = after_dilution(2:NE+2)/3.5e6;





tspan  = time_free_phages;
solved2 = ode45(@(t,y) one_step_eqn_before_dilution(t,y,theta,NE), tspan,after_dilution);
%y = solved2.y;
%time = solved2.x;

time = time_free_phages;
y = deval(solved2,time_free_phages);

end

