function [f,g] = VanderPolFun(t,x,mu)

f = [x(2); mu*(1-x(1)*x(1))*x(2)-x(1)];
g = x;