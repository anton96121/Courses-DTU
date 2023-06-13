function [f,J] = CSTR1D(t,x,param)
[F] = param(1);
Tin = param(2);
V = param(3);
beta = param(4);
Ca_in = param(5);
Cb_in = param(6);

Ca = Ca_in + (1/beta)*(Tin - x(1));
Cb = Cb_in + (2/beta)*(Tin - x(1));

k0 = exp(24.6);
EaR = 8500;

K_T = k0*exp(-EaR*(1/x(1)));
r = K_T*Ca*Cb; 
Rt = beta*r; 

f = zeros(1,1);
f(1) = (F/V)*(Tin - x(1)) + Rt;

% FIND JACOBIAN 
J = zeros(1,1);
a = x(1)^3+((-(Ca_in/2)-(Cb_in/4))*beta+(EaR/2)-Tin)*x(1)^2;
b = -(((Ca_in/2)+(Cb_in/4))*beta+Tin)*EaR*x(1);
c = +((EaR*(Ca_in*beta+Tin)*(((Cb_in*beta)/2)+Tin))/2);
J(1) = (1/V*beta*x(1)^2)*(4*V*k0*(a+b+c))*exp(-EaR*(1/x(1))-F*beta*x(1)^2);

end 