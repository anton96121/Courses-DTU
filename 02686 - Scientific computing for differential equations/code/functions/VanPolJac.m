function [J,dg] = VanPolJac(t,x,param)
%VANDERPOLF Jacobian of the Van der Pol problem
J = zeros(2,2);
J(1,1) = 0;
J(1,2) = 1;
J(2,1) = -2*param(1)*x(2)*x(1) - 1;
J(2,2) = param(1)*(1-x(1)^2);
dg = eye(2);
end

