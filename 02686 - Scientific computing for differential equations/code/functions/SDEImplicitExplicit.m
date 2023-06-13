function [x, t,function_calls,hs] = SDEImplicitExplicit(x0, f,jac, g, h, t0, T, param,paths, newtonTolerance, newtonMaxiterations,seed)
    N = ceil((T-t0)/h);
    hs = ones(1,N+1);
    hs = h.*hs;
    t = zeros(1,N+1);
    dim = length(x0);
    x = zeros(dim,N+1);
    t(1) = t0;
    x(:,1) = repmat(x0,1);
    function_calls = 0;
    
    W = brownian_motion(h,N+1,dim,seed);
    dt = f(t0, x(:,1), param);
    for k=2:N+1        
        dw = g(t(k-1), x(:,k-1), param);
        dW = W(:,k)-W(:,k-1);
        Rterm = x(:,k-1) + dw.*dW;
        t(k) = t(k-1)+h;
        xguess = Rterm + dt*h;
        [x(:,k), dt, function_calls] = newtonsMethod(f, jac,  t(:,k-1), Rterm, h, xguess, newtonTolerance, newtonMaxiterations, param, function_calls);
        function_calls = function_calls+2;
    end
    x = x';
    t = t';

end