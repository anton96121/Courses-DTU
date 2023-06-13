function [dx,x] = Lstable(t,x,param)
dx = -2000*(x-cos(t));
end

