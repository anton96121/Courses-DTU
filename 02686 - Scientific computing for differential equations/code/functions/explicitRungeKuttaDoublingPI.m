function [x,t,function_calls,hs,rs] = explicitRungeKuttaDoublingPI(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,initial_step_algo)
    N = ceil((T-t0)/hmin);
    s = length(b);
    d = length(x0);
    t = zeros(1,N+1);
    hs = zeros(2,N+1);
    rs = zeros(2,N+1);
    x = zeros(length(x0),N+1);
    k = zeros(d,s)';
    fac = 1;
    t(1) = t0;
    x(:,1) = x0;
    accept_step = false;
    function_calls = 0;
    p = p;
    
    if initial_step_algo % slide 4b, 9 whic is from p169 Hairer
        [h,function_calls] = initialStepSize(f,param,t0,x0);
    end
    
    i = 2;
    l = 2;
    runs = 1;
    while t(i-1)<=T
       while(~accept_step)
           if t(i-1)+h>T
               h = max(hmin,T-t(i-1));
           end
           
           k = 0*k;
           for j = 1:s
               k(j,:) = f(t(i-1)+h*c(j),x(:,i-1)+sum(h*k.*A(j,:)',1)', param);
               function_calls = function_calls +1;
           end
           oneAhead = x(:,i-1)+sum(h*k.*b,1)';

           hhalf = h/2;
           k = 0*k;
           for j = 1:s
               k(j,:) = f(t(i-1)+hhalf*c(j),x(:,i-1)+sum(hhalf*k.*A(j,:)',1)', param);
               function_calls = function_calls +1;
           end
           halfAhead = x(:,i-1)+sum(hhalf*k.*b,1)';

           k = 0*k;
           for j = 1:s
               k(j,:) = f(t(i-1)+2*hhalf*c(j),halfAhead+sum(hhalf*k.*A(j,:)',1)', param);
               function_calls = function_calls +1;
           end
           halfAhead = halfAhead+sum(hhalf*k.*b,1)';
           
           % Step doubling, see section II.4(s164) in Harier and slide 4c,3
           %sc = Atol + max(abs(halfAhead),abs(oneAhead))*Rtol;

           %e = norm((oneAhead-halfAhead)./sc,'Inf');
           e = abs(oneAhead-halfAhead);
           e = max(e./max(Atol,abs(halfAhead).*Rtol)); % r i john algo

           if e<=1
               x(:,i) = halfAhead;
               hs(1,i-1) = h;
               hs(2,i-1) = runs;
               rs(1,l-1) = e;
               rs(2,l-1) = 1;
               runs = 1;
               t(i) = t(i-1)+h;
               accept_step = true;
               if i==2
                   h = h*max(hmin,min(hmax,fac*(eps_tol/e)^(1/(p+1)))); % eq 4.13 in Harier  
                   e1 = max(e,10^(-10));
               else
                   %h = h*max(hmin,min(hmax,fac*(eps_tol/e)^(0.4/(p+1))*(e1/e)^(0.3/(p+1)))); % eq 4.13 in Harier
                   h = h*max(hmin,min(hmax,fac*(eps_tol/e)^(0.8/(p+1))*(e1/e)^(0.31/(p+1)))); % eq 4.13 in Harier
                   e1 = max(e,10^(-10));
               end
               
           else
               runs = runs +1;
               accept_step = false;
               rs(1,l-1) = e;
               rs(2,l-1) = 2;
               h = h*max(hmin,min(hmax,fac*(eps_tol/e)^(1/(p+1)))); % eq 4.13 in Harier  
               l = l+1;
           end

           
       end
       accept_step = false;
       i = i+1;
       l = l+1;
    end
    x = x(:,1:i-1)';
    t = t(1:i-1);
    hs = hs(1:2,1:i-2);
    rs = rs(1:2,1:l-2);
end
