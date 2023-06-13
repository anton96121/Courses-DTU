function [x,lambda]=rangeSpaceSolver(H,g,A,b)
Hg = H\g;
HA = H\A;
lambda = (A'*HA)\(b+A'*Hg);
x = HA*lambda-Hg;
end
