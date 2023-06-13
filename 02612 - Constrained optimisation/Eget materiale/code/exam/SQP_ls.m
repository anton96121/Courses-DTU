function [x,z,Hist] = SQP_ls(x0,obj,con,l,u,cl,cu,log, subsolver,precision, nonmonotone)
% SQP_ls       A sequential quadratic programing algorithm with a damped BFGS
%               update of the hessian and line search
%
%            min   f(x)
%             x
%            s.t   gu>= c(x) >= gl
%                   u>=  x >= l
%
%
% Syntax: [x,z,Hist] = SQP_ls(x0,obj,con,l,u,cl,cu,log, subsolver,precision, nonmonotone)
%
%         x             : Solution
%         z             : Lagrange multipliers
%         Hist          : Hist object with algorithm run-time information

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%  
    max_iter = 50;
    n = length(x0);
    epsilon = 10^(-precision);
    d = [l;-u;cl;-cu];
    mu = 0;
    functionCalls = 0;
    
  
    % Define relevant functions call
    x = x0;
    [f,df] = feval(obj,x);
    [c,dc] = feval(con,x);
    B = eye(n);
    m = size(c,1);
    z = ones(2*m+2*n,1);
    lid = 1:m;
    uid = (m+1):(2*m);
    clid = (2*m+1):(2*m+n);
    cuid = (2*m+n+1):(2*(n+m));
    
    % Define options for quadprog
    options = optimset('Display', 'off');
    
    % Log objects
    if log
        xHist = zeros(n,max_iter+1);
        pkHist = zeros(n,max_iter);
        timePerformence = zeros(1,max_iter);
        stepLength = zeros(1,max_iter);
        functionCalls = 2;
        
        xHist(:,1) = x0;
    end
    
    %Start for loop
    for i = 1:max_iter
        
        % Update lower and upper bounds for the quadrastart = cputime; approximation
        lk = -x+l;
        uk = -x+u;
        clk = -c+cl;
        cuk = -c+cu;
        
        % Solves local QP program
        if subsolver
            start = cputime;
            [pk,zhat] = intSQP(B,df,dc,lk,uk,clk,cuk,x);
            time = cputime-start;
        else
            start = cputime;
            [pk,~,~,~,lambda] = quadprog(B,df,-[dc'; -dc'],-[clk;-cuk],[],[],lk,uk,[], options);
            time = cputime-start;

            zhat = [lambda.lower; lambda.upper; lambda.ineqlin];
        end
            

        if any(isempty(pk) | isnan(pk) == true)
            disp(i)
            error('The program is infeasible. Try with infeasibility handling')
        end

        % Line search
        alpha = 1;
        pz = zhat-z;
        [c_l] = feval(con,x);
        functionCalls = functionCalls +1;
        c_ls = [x; -x; c_l; -c_l]-d;
        mu = max(abs(z),1/2*(mu+abs(z)));
        phi0 = phi(f,mu,c_ls);
        Dphi0 = dphi(df,pk,mu,c_ls); 
        
        while true
           x_ls = x + alpha*pk;
           f_ls = obj(x_ls);
           functionCalls = functionCalls +1;
           [c_l] = feval(con,x_ls);
           functionCalls = functionCalls +1;
           c_ls = [x; -x; c_l; -c_l]-d;

           phi1 = phi(f_ls,mu,c_ls);

           if phi1 <= phi0 + 0.1*alpha*Dphi0
               break
           else
               a = (phi1-(phi0+Dphi0*alpha))/alpha^2;
               alpha_min = -Dphi0/(2*a);

               alpha = min(0.9*alpha, max(alpha_min, 0.1*alpha));
           end
        end

        % non monotone strategy
        if (all(round(alpha*pk,precision)==zeros(n,1))) && nonmonotone
            alpha = 1;
        end


        % Update the current point
        z = z + alpha*pz;
        x = x + alpha*pk;
        mu = z;

        % For the quasi Newton update  
        dL = df - (z(lid)-z(uid)+dc*z(clid)-dc*z(cuid));

        % Update values for next iteration
        [f,df] = feval(obj,x);
        [c,dc] = feval(con,x);
        functionCalls = functionCalls +2;

        % Quasi newton update of the hessian
        dL2 = df - (z(lid)-z(uid)+dc*z(clid)-dc*z(cuid));

        p = alpha*pk;
        q = dL2-dL;
        theta = 1;
        Bp = (B*p);
        pBp = p'*Bp;
        
        if p'*q < 0.2*pBp
            theta = (0.8*pBp)/(pBp-p'*q);
        end
        r = theta*q+(1-theta)*(Bp);
        B = B + r*r'/(p'*r) - Bp*Bp'/pBp;
        
        % log information
        if log
            pkHist(:,i) = pk;
            xHist(:,i+1) = x;
            timePerformence(1,i) = time;
            stepLength(1,i) = sqrt(sum((alpha*pk).^2));
        end
        
        % Check for convergence
        if norm(dL2, 'inf')<epsilon
            if log
                pkHist = pkHist(:,1:i);
                xHist = xHist(:,1:i+1);
                timePerformence = timePerformence(:,1:i);
                stepLength = stepLength(1,1:i);
                
                Hist = struct('xHist', xHist, 'pkHist', pkHist, 'timePerformence', timePerformence, 'Iterations' ,i, 'stepLength', stepLength, 'functionCalls', functionCalls);
            else
                Hist = struct();
            end
            
            return
        end
   end

end
% Functions for the line search algorithm
function [val] = phi(f,mu,c)
val = f+mu'*abs(min(0,c));
end
function [val] = dphi(df,pk,mu,c)
val = df'*pk-mu'*abs(min(0,c));
end

function [val] = phialt(f,mu,c,z)
c = abs(min(0,c));

val = f-z'*c+mu'*c.^2;
end

function [val] = dphialt(df,pk,pz,mu,c,dc,z)
c = abs(min(0,c));

valx = df-dc*z+dc*(mu.*c);
valz = -c;

val = [valx; valz]'*[pk;pz];
end
