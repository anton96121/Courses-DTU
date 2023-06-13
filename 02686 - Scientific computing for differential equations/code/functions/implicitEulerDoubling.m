function [x,t,function_calls,hs,rs] = implicitEulerDoubling(f,jac,param,h,t0,T,x0,Atol,Rtol,hmin,hmax,eps_tol,initial_step_algo,newtonTolerance, newtonMaxiterations)
    N = ceil((T-t0)/hmin);
    p = 1;
    t = zeros(1,N+1);
    hs = zeros(2,N+1);
    rs = zeros(2,N+1);
    x = zeros(length(x0),N+1);
    t(1) = t0;
    x(:,1) = x0;
    accept_step = false;
    function_calls = 0;
    fac = 1;
    
    if initial_step_algo % slide 4b, 9 whic is from p169 Hairer
        [h,function_calls] = initialStepSize(f,param,t0,x0);
    end
    
    i = 2;
    l = 2;
    runs = 1;
    %h_old = h;
    xdot = f(t(1), x(:,1), param); % first guess
    function_calls = 3;
    while t(i-1)<=T
       while(~accept_step)
           if t(i-1)+h>T
               h = max(hmin,T-t(i-1));
           end
           
           t(i) = t(i-1)+h;
           xguess = x(:,i-1)+xdot*h; % We guess on a explicit euler step
           [oneAhead, ~, function_calls] = newtonsMethod(f, jac,  t(:,i-1), x(:,i-1), h, xguess, newtonTolerance, newtonMaxiterations, param, function_calls);

           hhalf = h/2;
           xguess = x(:,i-1)+xdot*hhalf; % We guess on a explicit euler step
           [halfAhead, xdotHalf, function_calls] = newtonsMethod(f, jac,  t(:,i-1), x(:,i-1), hhalf, xguess, newtonTolerance, newtonMaxiterations, param, function_calls);

           xguess = halfAhead+xdotHalf*hhalf; % We guess on a explicit euler step
           [halfAhead, xdot, function_calls] = newtonsMethod(f, jac,  t(:,i-1)+hhalf, halfAhead, hhalf, xguess, newtonTolerance, newtonMaxiterations, param, function_calls);
           
           % Step doubling, see section II.4(s164) in Harier and slide 4c,3
           e = abs(oneAhead-halfAhead);
           r = max(e./max(Atol,abs(halfAhead).*Rtol));
           
           if r<=1
               x(:,i) = halfAhead;
               hs(1,i-1) = h;
               hs(2,i-1) = runs;
               rs(1,l-1) = r;
               rs(2,l-1) = 1;
               runs = 1;
               t(i) = t(i-1)+h;
               accept_step = true;
               h = h*max(hmin,min(hmax,fac*(eps_tol/r)^(1/(1+p)))); % eq 4.13 in Harier
           else
               runs = runs +1;
               accept_step = false;
               rs(1,l-1) = r;
               rs(2,l-1) = 2;
               h = h*max(hmin,min(hmax,fac*(eps_tol/r)^(1/(1+p)))); % eq 4.13 in Harier
               l = l+1;
           end
           
       end
       accept_step = false;
       i = i+1;
       l = l+1;
    end
    x = x(:,1:i-1)';
    t = t(1:i-1);
    hs = hs(:,1:i-2);
    rs = rs(:,1:l-2);
end
