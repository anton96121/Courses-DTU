function [h,function_calls] = initialStepSize(f,param,t0,x0)
% slide 4b, slide 9 whic is from p169 Hairer
       % step a
       d0 = norm(x0,2);
       diff0 = f(t0,x0, param);
       d1 = norm(diff0,2);
       
       % step b
       if d0<= 10^(-5) || d1<= 10^(-5)
           h0 = 10^(-6);
       else
           h0 = 0.01*(d0/d1);
       end
       
       % step c
       x1 = x0+h0*diff0;
       diff1 = f(t0+h0,x1, param);
       
       % step d
       d2 = norm((diff1-diff0)./h0,2);
       
       % step e
       if max(d1,d2)<10^(-15)
           h1 = max(10^(-6),h0*10^(-3));
       else
           h1 = sqrt(0.01/max(d1,d2));
       end
       
       % step f
       h = min(100*h0,h1);
       function_calls = 2;
end
    