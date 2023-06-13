function [x,lambda, feasibility,working_set_save,x_save,p_save,v_save,r_save,lambda_save,step_save] = dualActiveSetMethod(H,g,A,b,max_iter,tol)
if nargin<7
  max_iter = 100;
end
if nargin<8
  tol = 0.00001;
end

x = -H\g;
[n,m] = size(A);

working_set_save = zeros(m,max_iter+1);
x_save = zeros(n,max_iter+1);
lambda_save = zeros(m,max_iter+1);
p_save = zeros(n,max_iter+1);
v_save = zeros(m,max_iter+1);
step_save = zeros(2,max_iter+1);
r_save =  zeros(1,max_iter+1);

working_set = zeros(m,1);
lambda = zeros(m,1);
Wk = find(working_set == 1);
WDk = find(working_set == 0);

x_save(:,1) = x;
working_set_save(:,1) = working_set;

for i=1:max_iter
    if ~isempty(Wk) && all(c(A,b,x,WDk)>0)
        working_set_save = working_set_save(:,1:i);
        x_save = x_save(:,1:i);
        p_save = p_save(:,1:i-1);
        v_save = v_save(:,1:i-1);
        r_save = r_save(:,1:i-1);
        lambda_save = lambda_save(:,1:i-1);
        step_save = step_save(:,1:i-1);
        
        
        feasibility = 1;
        return
    end
    cD = c(A,b,x,WDk);
    r = cD==max(cD(cD<0));
    r_true = WDk(r);
    
    
    while c(A,b,x,r_true)<-tol
        r_save(:,i) = r_true;

        [p,v] = rangeSpaceSolver(H,-A(:,r_true),A(:,Wk),zeros(length(Wk),1));
        
        p_save(:,i) = p;
        v_save(Wk,i) = v;
        
        if norm(A(:,r_true)'*p)<tol
            if all(v>=-tol)
                
                feasibility = 0;
                return
            else
                smaller = find(v < 0);
                t = -lambda(Wk(smaller))./v(smaller);
                j = min(t)==t;
                j = smaller(j);
                j = Wk(j);
                t = min(t);
                
                lambda(Wk) = lambda(Wk) + t*v;
                lambda(r_true) = lambda(r_true) + t;
                working_set(j) = 0;
                Wk = find(working_set == 1);
                WDk = find(working_set == 0);
                
                
                working_set_save(:,i+1) = working_set;
                x_save(:,i+1) = x;
                lambda_save(:,i) = lambda;
                step_save(2,i) = t;
            end
        else
            if all(v>=0)
                tmax = Inf;
            else
                smaller = find(v < 0);
                tmax = -lambda(Wk(smaller))./v(smaller);
                j = min(tmax)==tmax;
                j = smaller(j);
                j = Wk(j);
                tmax = min(tmax);
            end
            tstar = -c(A,b,x,r_true)/(A(:,r_true)'*p);
            if tstar <= tmax
                x = x+tstar*p;
                lambda(Wk) = lambda(Wk) + tstar*v;
                lambda(r_true) = lambda(r_true)+tstar;
                working_set(r_true) = 1;
                Wk = find(working_set == 1);
                WDk = find(working_set == 0);
                
                working_set_save(:,i+1) = working_set;
                x_save(:,i+1) = x;
                lambda_save(:,i) = lambda;
                step_save(1,i) = tstar;
                step_save(2,i) = tmax;
            else
                x = x + tmax*p;
                lambda(Wk) = lambda(Wk)+tmax*v;
                lambda(r_true) = lambda(r_true)+tmax;
                working_set(j) = 0;
                Wk = find(working_set == 1);
                WDk = find(working_set == 0);
                
                working_set_save(:,i+1) = working_set;
                x_save(:,i+1) = x;
                lambda_save(:,i) = lambda;
                step_save(1,i) = tstar;
                step_save(2,i) = tmax;
            end
        end
    end
end
end

