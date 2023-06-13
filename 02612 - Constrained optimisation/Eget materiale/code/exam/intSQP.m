function [x,z,feasible,i] = intSQP(B,df,dc,lk,uk,clk,cuk,x0)
% intSQP   An interior point solver based on Mehrota's predictor-corrector
%           primal-dual interior point algorithm. It takes 
%            problems of the form
%
%            min   x'*H*x+g'x
%             x
%            s.t   gu>= cx >= gl
%                   u>=  x >= l
%
%
% Syntax: [x,z,feasible,i] = intSQP(B,df,dc,lk,uk,clk,cuk,x0)
%
%         x             : Solution
%         z             : Lagrange multipliers
%         feasible      : Flag to indicate feasibility
%         i             : Iterations used

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    % Sets constants for the algorithm
    warning('off','all')    
    n = length(x0);
    m = length(cuk);
    mc = n*2+2*m;
    epsilon = 0.0001;
    max_iter = 100;
    eta = 0.995;
    feasible = 1;
    
    % Initial values
    x = x0;
    z = ones(mc,1);
    s = ones(mc,1);
    d = [lk;-uk;clk;-cuk];
    e = ones(2*n+2*m,1);
    
    % Makes sure non of the following matrix operations are singular
    while(any(s==0))
        x = x+0.000001;
        s = [x;-x; dc'*x; -dc'*x]-d;
    end
    
    %Initilize constraint specific slacks and lagrange multipliers
    sl = s(1:n);
    su = s(n+1:n*2);
    scl = s(2*n+1:2*n+m);
    scu = s(m+2*n+1:n*2+2*m);
    zl = z(1:n);
    zu = z(n+1:n*2);
    zcl = z(2*n+1:2*n+m);
    zcu = z(m+2*n+1:n*2+2*m);
    
    %initial residuals
    rL = B*x+df-(zl-zu+dc*zcl-dc*zcu);
    rC = s+d - [x; -x; dc'*x; -dc'*x];
    
    % Start point heuristic
    zsl = diag(zl./sl);
    zsu = diag(zu./su);
    zslc = zcl./scl;
    zsuc = zcu./scu;
    zc = zslc + zsuc;
    Hbar = B + zsl + zsu + bsxfun(@times,zc',dc)*dc';
    [L,D,P] = ldl(Hbar,'lower');

    % Affine step for the start point heuristic
    rCs = (rC-s);
    rLbar = rL - [ zsl -zsu bsxfun(@times,zslc',dc) bsxfun(@times,zsuc',-dc)]*rCs;
    
    rhs = -rLbar;
    dxAff = P*(L' \ (D \ (L \ (P'*rhs) )));

    dzAff = - [ zsl; -zsu; bsxfun(@times,zslc,dc'); bsxfun(@times,zsuc,-dc')]*dxAff + (z./s).*rCs;
    dsAff = -s-(s./z).*dzAff;
    
    %Update of starting point
    z = max(1,abs(z+dzAff));
    s = max(1,abs(s+dsAff));
    
    %Update of initial residuals
    sl = s(1:n);
    su = s(n+1:n*2);
    scl = s(2*n+1:2*n+m);
    scu = s(m+2*n+1:n*2+2*m);
    zl = z(1:n);
    zu = z(n+1:n*2);
    zcl = z(2*n+1:2*n+m);
    zcu = z(m+2*n+1:n*2+2*m);

    
    rL = B*x+df-(zl-zu+dc*zcl-dc*zcu);
    rC = s+d - [x; -x; dc'*x; -dc'*x];
    rSZ = s.*z;
    
    % Initial dual gap
    dualGap = (z'*s)/(mc);
    dualGap0 = dualGap;
    
    for i = 1:max_iter
        zsl = diag(zl./sl);
        zsu = diag(zu./su);
        zslc = zcl./scl;
        zsuc = zcu./scu;
        zc = zslc + zsuc;
        
        % Factorization
        Hbar = B + zsl + zsu + bsxfun(@times,zc',dc)*dc';
        [L,D,P] = ldl(Hbar,'lower');

        % Affine step
        rCs = (rC-s);
        rLbar = rL - [ zsl -zsu bsxfun(@times,zslc',dc) bsxfun(@times,zsuc',-dc)]*rCs;

        rhs = -rLbar;
        dxAff = P*(L' \ (D \ (L \ (P'*rhs) )));

        dzAff = - [ zsl; -zsu; bsxfun(@times,zslc,dc'); bsxfun(@times,zsuc,-dc')]*dxAff + (z./s).*rCs;
        dsAff = -s-(s./z).*dzAff;
        
        % Compute max alpha affine
        dZS = [dzAff; dsAff];
        alphas = (-[z;s]./dZS);
        alphaAff = min([1;alphas(dZS<0)]);
           
        dualGapAff = ((z+alphaAff*dzAff)'*(s+alphaAff*dsAff))/(mc);
        sigma = (dualGapAff/dualGap)^3;
        
        
        % Affine-Centering-Correction Direction
        rSZbar = rSZ + dsAff.*dzAff-sigma*dualGap*sigma*e;        
        rLbar = rL - [ zsl -zsu bsxfun(@times,zslc',dc) bsxfun(@times,zsuc',-dc)]*(rC-rSZbar./z);
        
        rhs = -rLbar;
        dx = P*(L' \ (D \ (L \ (P'*rhs) )));
        
        dz = - [ zsl; -zsu; bsxfun(@times,zslc,dc'); bsxfun(@times,zsuc,-dc')]*dx + (z./s).*(rC-rSZbar./z);
        ds = -rSZbar./z-(s./z).*dz;
        
        %compute max alpha
        dZS = [dz; ds];
        alphas = (-[z;s]./dZS);
        alpha = min([1;alphas(dZS<0)]);
        
        alphaBar = eta*alpha;
        
        % Update of position 
        x = x + alphaBar * dx;
        z = z + alphaBar * dz;
        s = s + alphaBar * ds;
        
        % Update of residuals
        sl = s(1:n);
        su = s(n+1:n*2);
        scl = s(2*n+1:2*n+m);
        scu = s(m+2*n+1:n*2+2*m);
        zl = z(1:n);
        zu = z(n+1:n*2);
        zcl = z(2*n+1:2*n+m);
        zcu = z(m+2*n+1:n*2+2*m);


        rL = B*x+df-(zl-zu+dc*zcl-dc*zcu);
        rC = s+d - [x; -x; dc'*x; -dc'*x];
        rSZ = s.*z;
        
        % Compute the dual gap
        dualGap = (z'*s)/(mc);
        
        if any(isnan(z))
            feasible = 0;
            return
        end
        
        % Check for convergence
        if(dualGap <= epsilon*0.01*dualGap0)
           return 
        end
    end  
end

