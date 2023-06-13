function [dx,x] = Astable_half(t,x,param)
dx = zeros(2,1);
dx(1) = x(2);
dx(2) = -x(1);
end

