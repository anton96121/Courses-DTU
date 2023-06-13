function [x,t,function_calls,hs] = explicitRungeKutta(f,param,h,t0,T,x0,A,b,c)
    N = ceil((T-t0)/h);
    t = zeros(1,N+1);
    hs = ones(1,N+1);
    hs = h.*hs;
    x = zeros(length(x0),N+1);
    s = length(b);
    k = zeros(length(x0),s)';
    t(1) = t0;
    x(:,1) = x0;
    function_calls = 0;
    
    Ah = h*A;
    bh = h*b;
    ch = h*c;
    
    for i = 2:N+1
       if t(i-1)+h > T
            h = T-t(i-1);
       end
       k = 0*k;
       for j = 1:s
           k(j,:) = f(t(i-1)+ch(j),x(:,i-1)+sum(k.*Ah(j,:)',1)', param);
           function_calls = function_calls +1;
       end
       x(:,i) = x(:,i-1)+sum(k.*bh,1)';
       t(i) = t(i-1)+h;
    end
    x = x';
end