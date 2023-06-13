function [x,z,Hist] = SQP_trust(x0,obj,cons,l,u,cl,cu,log,precision,trust_region, penalty)
% SQP_trust       A sequential quadratic programing algorithm with a damped BFGS
%                   update of the hessian, line search and infeasibility handling
%
%            min   f(x)
%             x
%            s.t   gu>= c(x) >= gl
%                   u>=  x >= l
%
%
% Syntax: [x,z,Hist] = SQP_trust(x0,obj,cons,l,u,cl,cu,log,precision,trust_region, penalty)
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
    mu = penalty;
    dk0 = trust_region;
    dk = dk0;
    time = 0;
    dL2 = Inf;

    % Define relevant functions call
    x = x0;
    [f,df] = obj(x);
    [c,dc] = cons(x);
    B = eye(n);
    m = size(c,1);
    mu_vec = mu*ones(n*2+m*2,1);
    z = ones(2*m+2*n,1);
    lid = 1:m;
    uid = (m+1):(2*m);
    clid = (2*m+1):(2*m+n);
    cuid = (2*m+n+1):(2*(n+m));
    
    % Log objects
    if log
        xHist = zeros(n,max_iter+1);
        pkHist = zeros(n,max_iter);
        timePerformence = zeros(1,max_iter);
        trustRegion = zeros(1,max_iter+1);
        rhos = zeros(1,max_iter);
        functionCalls = 2;
        
        xHist(:,1) = x0;
        trustRegion(:,1) = dk;
    end
    
    % Define options for quadprog
    options = optimset('Display', 'off');
    

    % Start for loop
    for i = 1:max_iter
        
        % Update lower and upper bounds for the quadrastart = cputime; approximation
        lk = -x+l;
        uk = -x+u;
        clk = -c+cl;
        cuk = -c+cu;
        
        % Define infesibility program
        Hinf = zeros(3*n+2*m);
        Hinf(1:n,1:n) = B;
        Cinf = [ eye(n) -eye(n) dc -dc zeros(n,2*m+2*n) eye(n)  -eye(n);eye(2*n+2*m) eye(2*n+2*m) zeros(2*m+2*n,n*2)]';
        ginf = [df; mu_vec];
        dinf = [lk; -uk;clk; -cuk; zeros(2*m+2*n,1); -dk*ones(2*m,1)];
        start = cputime;
        [pk,~,~,~,zhat] = quadprog(Hinf,ginf,-Cinf,-dinf,[],[],[],[],[],options);
        time = time + cputime-start;

        zhat = zhat.ineqlin(1:2*m+2*n);
        pk = pk(1:n);
        
        %Update penalty
        zinf = vecnorm(z,'Inf');
        mu = max(1/2*(mu+zinf),zinf);
        mu_vec = mu*ones(n*2+m*2,1);
        
        %Calculate relevant values for the trust region
        [c_full, dc_full] = cons(x,true,dinf(1:2*n+2*m));
        functionCalls = functionCalls +1;

        [c_p_full, ~] = cons(x+pk,true,dinf(1:2*n+2*m));
        functionCalls = functionCalls +1;
        
        qp0 = f+df'*pk+1/2*pk'*B*pk+mu_vec'*max(0,-(c_full+dc_full'*pk));
        q0 = f+mu_vec'*max(0,-(c_full));
        phi1 = q0;
        
        [f_p,~] = obj(x+pk);
        functionCalls = functionCalls +1;
        phi1p = f_p+mu_vec'*max(0,-(c_p_full));
        
        rho = (phi1-phi1p)/(q0-qp0);
        
        gamma = min(max((2*rho-1)^3+1,0.25),2);
        
        % If trust region is accepted
        if rho>0
            % Update the current point
            z = zhat;
            x = x+pk;
        

            % For the quasi Newton update  
            dL = df - (z(lid)-z(uid)+dc*z(clid)-dc*z(cuid));



            % Update values for next iteration
            [f,df] = feval(obj,x);
            [c,dc] = feval(cons,x);

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
            
            % Update the trust region
            dk = gamma*dk;
        
        % If the trust region is not accepted
        else
            % Update the trust region
            dk = gamma * vecnorm(pk,'Inf');
        end

        % Store relevant logging information
        if log
            pkHist(:,i) = pk;
            xHist(:,i+1) = x;
            timePerformence(1,i) = time;
            trustRegion(:,i+1) = dk;
            rhos(:,i) = rho;
            time = 0;
        end

        % Check for convergence 
        if norm(dL2, 'inf')<epsilon
            if log
                pkHist = pkHist(:,1:i);
                xHist = xHist(:,1:i+1);
                timePerformence = timePerformence(:,1:i);
                rhos = rhos(:,1:i);
                trustRegion = trustRegion(:,1:i+1);
                
                
                Hist = struct('xHist', xHist, 'pkHist', pkHist, 'timePerformence', timePerformence, 'Iterations' ,i, 'functionCalls', functionCalls, 'rho', rhos, 'trustRegion',trustRegion);
            else
                Hist = struct();
            end
            
            return
        end
   end

end

