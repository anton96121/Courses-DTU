function [x,t,function_calls,hs,rs] = explicitRungeKuttaEmbeddedPI(f,param,h,t0,T,x0,A,b,bhat,c,Atol,Rtol,hmin,hmax,eps_tol,p,initial_step_algo)
    N = ceil((T-t0)/hmin);
    s = length(b);
    d = length(x0);
    t = zeros(1,N+1);
    hs = zeros(2,N+1);
    rs = zeros(2,N+1);
    fac = 1;
    x = zeros(length(x0),N+1);
    k = zeros(d,s)';
    t(1) = t0;
    x(:,1) = x0;
    accept_step = false;
    function_calls = 0;
    
    
    if initial_step_algo % slide 4b, 9 whic is from p169 Hairer
        [h,function_calls] = initialStepSize(f,param,t0,x0);
    end
    
    runs = 1;
    i = 2;
    l = 2;
    while t(i-1)<=T
        
       while(~accept_step)
           k = 0*k;
           for j = 1:s
               k(j,:) = f(t(i-1)+h*c(j),x(:,i-1)+sum(h*k.*A(j,:)',1)', param);
               function_calls = function_calls +1;
           end
           oneAhead = x(:,i-1)+sum(h*k.*b,1)';
           oneAheadhat = x(:,i-1)+sum(h*k.*bhat,1)';
           
           % Step doubling, see section II.4(s164) in Harier and slide 4c,3
           e = abs(oneAhead-oneAheadhat);
           e = max(e./max(Atol,abs(oneAhead).*Rtol)); % r i john algo

           if e<=1
               x(:,i) = oneAhead;
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
               rs(1,l-1) = e;
               rs(2,l-1) = 2;
               accept_step = false;
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
    hs = hs(:,1:i-2);
    rs = rs(:,1:l-2);
end