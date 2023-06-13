function [dx,g] = VanPol(t,x,param)
dx = zeros(2,1);
dx(1) = x(2);
dx(2) = param(1)*(1-x(1)^2)*x(2)-x(1);

dx  = repmat(dx,length(x)/2,1);
g = x;
end
