function [x,lambda,feasibility,working_set_save,x_save,p_save,lambda_save,step_save] = primalActiveSetMethod(H,g,A,b,x0,W0,max_iter,tol)
if nargin<7
  max_iter = 100;
end
if nargin<8
  tol = 0.00001;
end

[n,m] = size(A);
working_set_save = zeros(m,max_iter+1);
x_save = zeros(n,max_iter+1);
lambda_save = zeros(m,max_iter+1);
p_save = zeros(n,max_iter+1);
step_save = zeros(1,max_iter+1);

lambdas = zeros(m,1);
working_set = zeros(m,1);
working_set(W0) = 1;
x = x0;
x_save(:,1) = x0;
working_set_save(:,1) = working_set;
feasibility = 0;

for i=1:max_iter
   Wk = find(working_set == 1);
   %Solve direction
   [p,lambda] = rangeSpaceSolver(H,(H*x+g),A(:,Wk),zeros(length(Wk),1));
   
   p_save(:,i) = p;
   lambda_save(Wk,i) = lambda;
   if norm(p)<tol
       if ~any(lambda<=0)
           lambdas(Wk) = lambda;
           lambda = lambdas;
           x_save(:,i+1) = x;
           working_set_save(:,i+1) = working_set;
           
           x_save = x_save(:,1:i+1);
           p_save = p_save(:,1:i);
           lambda_save = lambda_save(:,1:i);
           step_save = step_save(:,1:i);
           working_set_save = working_set_save(:,1:i+1);
           feasibility = 1;
           return
       else
           j = lambda==min(lambda(lambda<0));
           j = Wk(j);
           working_set(j) = 0;
           
           x_save(:,i+1) = x;
           working_set_save(:,i+1) = working_set;
       end
   else
       notWk = find(working_set == 0);
       denom = (A(:,notWk)'*p);
       dummy = (b(notWk(denom<-tol))-A(:,notWk(denom<-tol))'*x)./denom(denom<-tol);
       alpha = min(1,min(dummy));
       
       f = (find(min(dummy)==dummy));
       k=find((denom<0));
       j = k(f(1));
       
       step_save(:,i) = alpha;
       
       if alpha ~= 1
           x = x+alpha*p;
           working_set(j) = 1;
           x_save(:,i+1) = x;
           working_set_save(:,i+1) = working_set;
       else
           x = x+p;
           x_save(:,i+1) = x;
           working_set_save(:,i+1) = working_set;
       end
   end
end     
end

