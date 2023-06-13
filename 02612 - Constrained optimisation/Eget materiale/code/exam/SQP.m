function [x,z,Hist] = SQP(x0,obj,con,l,u,cl,cu,log, subsolver,precision)
% SQP       A sequential quadratic programing algorithm with a damped BFGS
%           update of the hessian
%
%            min   f(x)
%             x
%            s.t   gu>= c(x) >= gl
%                   u>=  x >= l
%
%
% Syntax: [x,z,Hist] = SQP(x0,obj,con,l,u,cl,cu,log, subsolver,precision)
%
%         x             : Solution
%         z             : Lagrange multipliers
%         Hist          : Hist object with algorithm run-time information

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    max_iter = 100;
    n = length(x0);
    epsilon = 10^(-precision);
    functionCalls = 0;

    % Define relevant functions call
    x = x0;
    [~,df] = feval(obj,x);
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
        functionCalls = 2;
        
        xHist(:,1) = x0;
    end
    

    % Start for loop
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
            error('The program is infeasible. Try with infeasibility handling')
        end

        pz = zhat-z;


        % Update the current point
        z = z + pz;
        x = x + pk;
        

        % For the quasi Newton update  
        dL = df - (z(lid)-z(uid)+dc*z(clid)-dc*z(cuid));

        % Update values for next iteration
        [~,df] = feval(obj,x);
        [c,dc] = feval(con,x);
        functionCalls = functionCalls +2;

        % Quasi newton update of the hessian
        dL2 = df - (z(lid)-z(uid)+dc*z(clid)-dc*z(cuid));

        p = pk;
        q = dL2-dL;
        theta = 1;
        Bp = (B*p);
        pBp = p'*Bp;
        
        if p'*q < 0.2*pBp
            theta = (0.8*pBp)/(pBp-p'*q);
        end
        r = theta*q+(1-theta)*(Bp);
        B = B + r*r'/(p'*r) - Bp*Bp'/pBp;
        
        if log
            pkHist(:,i) = pk;
            xHist(:,i+1) = x0;
        end
        
        % log information
        if log
            pkHist(:,i) = pk;
            xHist(:,i+1) = x;
            timePerformence(1,i) = time;
        end
        
        % Check for convergence
        if norm(dL2, 'inf')<epsilon
            if log
                pkHist = pkHist(:,1:i);
                xHist = xHist(:,1:i+1);
                timePerformence = timePerformence(:,1:i);
                
                Hist = struct('xHist', xHist, 'pkHist', pkHist, 'timePerformence', timePerformence, 'Iterations' ,i, 'functionCalls', functionCalls);
            else
                Hist = struct();
            end
            
            return
        end
   end

end

