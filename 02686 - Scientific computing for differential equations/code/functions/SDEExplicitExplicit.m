function [x, t,function_calls,hs] = SDEExplicitExplicit(x0, f, g, h, t0, T, params,paths,seed)
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
    
    for k=2:N+1
        dW = W(:,k)-W(:,k-1);
        dt = f(t(k-1), x(:,k-1), params);
        dw = g(t(k-1), x(:,k-1), params);
        x(:,k) = x(:,k-1) + dt*h + dw.*dW;
        t(k) = t(k-1)+h;
        function_calls = function_calls+2;
    end
    x = x';
    t = t';

end