function [dx,x_out] = Astable_half_Jac(t,x,param)
dx = [0 1; -1 0];
x_out = eye(2);
end