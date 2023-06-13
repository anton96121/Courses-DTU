function [dfdx,dgdx]=VanderPolJac(t,x,mu)

dfdx = [0 1;-2*mu*x(1)*x(2)-1 mu*(1-x(1)*x(1))];
dgdx = eye(2);