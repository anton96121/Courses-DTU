function [x,y,z,s, iter, ldltime] = primalDualInteriorMethod_box(H,g,A,b,l,u,x0,y0,z0,s0)
% primalDualInteriorMethod_box   An interior point solver based on Mehrota's predictor-corrector
%                                   primal-dual interior point algorithm. It takes 
%                                   problems of the form
%
%            min    x'Hx+g'x
%             x
%            s.t     Ax  = b
%                u>=  x >= l
%
%
% Syntax: [x,y,z,s, iter, ldltime] = primalDualInteriorMethod_box(H,g,A,b,l,u,x0,y0,z0,s0)
%
%         x             : Solution
%         y             : Equality lagrange multipliers
%         z             : Inequality lagrange multipliers
%         s             : Slack variables
%         iter          : Iterations used
%         ldltime       : Time used on ldl factorization

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    % Sets constants for the algorithm
    mIn = length(x0);
    m = length(y0);
    epsilon = 0.000001;
    ldltime = 0;
    max_iter = 100;
    eta = 0.995;
    iter = 0;
    
    % Initial values  
    x = x0;
    y = y0;
    z = z0;
    s = s0;
    
    % Makes sure non of the following matrix operations are singular
    while(any(s==0))
        x = x+0.000001;
        s = [x-l;-x+u];
    end
    
    %Initilize constraint specific slacks and lagrange multipliers 
    e = ones(mIn*2,1);
    sl = s(1:mIn);
    su = s(mIn+1:mIn*2);
    zl = z(1:mIn);
    zu = z(mIn+1:mIn*2);
    
    %initial residuals
    rL = H*x+g-A*y-(z(1:mIn)-z(mIn+1:mIn*2));
    rA = b-A'*x;
    rC = s+[l; -u] - [x; -x];
    
    % Start point heuristic
    zsl = zl./sl;
    zsu = zu./su;
    Hbar = H + diag(zsl + zsu);
    KKT = [Hbar -A; -A' zeros(m)];
    KKT = sparse(KKT);
    start = cputime;
    [L,D,p] = ldl(KKT,'vector');
    ldltime = ldltime + cputime-start;
    

    % Affine step
    rCs = (rC-s);
    rLbar = rL - zsl.*rCs(1:mIn) +zsu.*rCs(1+mIn:2*mIn);
    rhs = -[rLbar ; rA];
    
    solution(p) = L' \ (D \ (L \rhs(p)));

    dxAff = solution(1:length(x))';

    dzAff = - [zsl.*dxAff; -zsu.*dxAff] + (z./s).*rCs;
    dsAff = -s-(s./z).*dzAff;
    
    %Update of starting point
    z = max(1,abs(z+dzAff));
    s = max(1,abs(s+dsAff));
    
    %Update of initial residuals
    sl = s(1:mIn);
    su = s(mIn+1:mIn*2);
    zl = z(1:mIn);
    zu = z(mIn+1:mIn*2);

    
    rL = H*x+g-A*y-(zl-zu);
    rA = b-A'*x;
    rC = s+[l; -u] - [x; -x];
    
    % Initial dual gap
    dualGap = (z'*s)/(2*mIn);
    dualGap0 = dualGap;

    
    for i = 1:max_iter
        iter = iter + 1;
        zsl = zl./sl;
        zsu = zu./su;
        Hbar = H + diag(zsl + zsu);
        KKT = [Hbar -A; -A' zeros(m)];
        KKT = sparse(KKT);
        start = cputime;
        [L,D,p] = ldl(KKT,'vector');
        ldltime = ldltime + cputime-start;
        
        % Affine step
        rCs = (rC-s);
        rLbar = rL - zsl.*rCs(1:mIn) +zsu.*rCs(1+mIn:2*mIn);

        rhs = -[rLbar ; rA];
        solution(p) = L' \ (D \ (L \rhs(p)));
        
        dxAff = solution(1:length(x))';
        
        dzAff = - [zsl.*dxAff; -zsu.*dxAff] + (z./s).*rCs;
        dsAff = -s-(s./z).*dzAff;
        
        %compute max alpha affine
        dZS = [dzAff; dsAff];
        alphas = (-[z;s]./dZS);
        alphaAff = min([1;alphas(dZS<0)]);
           
        dualGapAff = ((z+alphaAff*dzAff)'*(s+alphaAff*dsAff))/(2*mIn);
        sigma = (dualGapAff/dualGap)^3;
        
        
        % Affine-Centering-Correction Direction
        rSZz =  s + dsAff.*dzAff./z-dualGap*sigma*e./z;
        rCs = (rC-rSZz);
        rLbar = rL - zsl.*rCs(1:mIn) +zsu.*rCs(1+mIn:2*mIn);
        
        rhs = -[rLbar ; rA];
        solution(p) = L' \ (D \ (L \rhs(p)));
        
        dx = solution(1:length(x))';
        dy = solution(length(x)+1:length(x)+length(y))';
        
        dz = - [zsl.*dx; -zsu.*dx] + (z./s).*rCs;
        ds = -rSZz-(s./z).*dz;
        
        %compute max alpha
        dZS = [dz; ds];
        alphas = (-[z;s]./dZS);
        alpha = min([1;alphas(dZS<0)]);
        
        alphaBar = eta*alpha;
        
        % Update of position 
        x = x + alphaBar * dx;
        y = y + alphaBar * dy;
        z = z + alphaBar * dz;
        s = s + alphaBar * ds;
        
        % Update of residuals 
        sl = s(1:mIn);
        su = s(mIn+1:mIn*2);
        zl = z(1:mIn);
        zu = z(mIn+1:mIn*2);


        rL = H*x+g-A*y-(z(1:mIn)-z(mIn+1:mIn*2));
        rA = b-A'*x;
        rC = s+[l; -u] - [x; -x];

        % Compute the dual gap        
        dualGap = (z'*s)/(2*mIn);
        
        % Check for convergence
        if(dualGap <= epsilon*0.01*dualGap0)
           return
        end
    end  
end

