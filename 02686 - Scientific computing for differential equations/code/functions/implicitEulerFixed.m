function [x,t,function_calls,hs] = implicitEulerFixed(f,jac,param,h,t0,T,x0,newtonTolerance,newtonMaxiterations)

    N = ceil((T-t0)/h);
    t = zeros(1,N+1);
    hs = ones(1,N+1);
    hs = h.*hs;
    x = zeros(length(x0),N+1);
    t(1) = t0;
    x(:,1) = x0;
    function_calls = 0;


% Fix 2:N+1
xdot = f(t(1), x(:,1), param); % first guess
function_calls = function_calls +1;
for i = 2:N+1
    if t(i-1)+h > T
        h = T-t(i-1);
    end
    t(i) = t(i-1)+h;
    xguess = x(:,i-1)+xdot*h; % We guess on a explicit euler step
    [x(:,i), xdot, function_calls] = newtonsMethod(f, jac,  t(:,i-1), x(:,i-1), h, xguess, newtonTolerance, newtonMaxiterations, param, function_calls);
end

t = t';
x = x';
end

