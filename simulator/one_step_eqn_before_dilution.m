function dydt = one_step_eqn_before_dilution(t,y,theta,NE)
%ONE_STEP_EQN Summary of this function goes here
%   Detailed explanation: Differential equations for SEIV (1,1,NE)


S = y(1);
E_mat(1:NE) = y(2:NE+1);
I = y(NE+2);
V = y(NE+3);

r = theta(1);
phi = theta(2);
tau = theta(3);
beta = theta(4);

etaeff = ((NE+1)/tau);

dotS = r*S - phi*V*S;
dotE1 = phi*S*V - etaeff * E_mat(1);
dotE_mat = etaeff.*E_mat(1:end-1) - etaeff.*E_mat(2:end);
dotI = etaeff * (E_mat(end) - I);
dotV = beta * etaeff * I - V*phi*(S+I+sum(E_mat));

dydt = [dotS,dotE1,dotE_mat,dotI,dotV]';

end

