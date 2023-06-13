function [f,J] = CSTR3D(t,x,param)
F = param(1);
Tin = param(2);
V = param(3);
beta = param(4);
Ca_in = param(5);
Cb_in = param(6);

EaR = 8500;
expa = exp(24.6); 

K_T = exp(24.6)*exp(-EaR*(1/x(3)));
r = K_T*x(1)*x(2); 
Ra = -r;
Rb = -2*r;
Rt = beta*r; 

f = zeros(3,1);
f(1) = (F/V)*(Ca_in - x(1)) + Ra;
f(2) = (F/V)*(Cb_in - x(2)) + Rb;
f(3) = (F/V)*(Tin - x(3)) + Rt;

% FIND JACOBIAN 
J = zeros(3,3);

J(1,1) = -(F/V)-expa*exp(-8500/x(3))*x(2);
J(1,2) = -expa*exp(-8500/x(3))*x(1);
J(1,3) = -((8500*expa*exp(-8500/x(3))*x(1)*x(2))/x(3)^2);

J(2,1) = -2*expa*exp(-8500/x(3))*x(2);
J(2,2) = -(F/V)-2*expa*exp(-8500/x(3))*x(1);
J(2,3) = -((1700*expa*exp(-8500/x(3))*x(1)*x(2))/x(3)^2);

J(3,1) = beta*expa*exp(-8500/x(3))*x(2);
J(3,2) = beta*exp(-8500/x(3))*x(1);
J(3,3) = -(F/V)+((8500*beta*expa*exp(-8500/x(3))*x(1)*x(2))/x(3)^2);

end 