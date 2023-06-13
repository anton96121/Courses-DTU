function [t,y] = itrapez(f,tspan,y0,N,fjac,atol,maxiter)

h = (tspan(2)-tspan(1))/N;
t = tspan(1)*ones(1,N+1) + h*(0:N);
d = length(y0);

t(1) = tspan(1);
y(:,1) = y0;
for i=2:N+1
yn = y(:,i-1); tn = t(i-1);
F = @(Z) Z - h*[zeros(d,1); f(tn,yn + Z(1:d))/2 + f(tn + h,yn + Z(d+1:2*d))/2];
dF = @(Z) eye(2*d) - h*[zeros(d,2*d);...
fjac(tn,yn + Z(1:d))/2, fjac(tn + h,yn + Z(d+1:2*d))/2];

Z = newton_itra(F,dF,zeros(2*d,1),atol,maxiter);

y(:,i) = yn + h*( f(tn,yn + Z(1:d))/2 + f(tn + h,yn + Z(d+1:2*d))/2 );
end
end

