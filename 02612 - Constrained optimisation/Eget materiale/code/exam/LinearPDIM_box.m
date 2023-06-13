function [x,y,z,s, iter] = LinearPDIM_box(g,A,b,l,u,x0,y0,z0,s0)
% LinearPDIM_box   An interior point solver based on Mehrota's predictor-corrector
%                   primal-dual interior point algorithm. It takes 
%                   problems of the form
%
%            min    g'x
%             x
%            s.t     Ax  = b
%                u>=  x >= l
%
%
% Syntax: [x,y,z,s, iter] = LinearPDIM_box(g,A,b,l,u,x0,y0,z0,s0)
%
%         x             : Solution
%         y             : Equality lagrange multipliers
%         z             : Inequality lagrange multipliers
%         s             : Slack variables
%         iter          : Iterations used

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    % Sets constants for the algorithm
    mIn = length(u);
    epsilon = 0.000000001;
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
    rL = g-A*y-(z(1:mIn)-z(mIn+1:mIn*2));
    rA = b-A'*x;
    rC = s+[l; -u] - [x; -x];
    
    % Start point heuristic
    zsl = zl./sl;
    zsu = zu./su;
    
    % Affine step for the start point heuristic
    rCs = (rC-s);
    rLbar = rL - zsl.*rCs(1:mIn) +zsu.*rCs(1+mIn:2*mIn);
    Hbar_diagonal_inverse = 1./(zsl+zsu);
    
    % Calculate the factor in the normal equation
    normalfactor = A' * (Hbar_diagonal_inverse .*  A);
    R = chol(normalfactor);
    
    mu_rhs = rA + A' * (Hbar_diagonal_inverse .* rLbar);
    dyAff = R \ (R' \ mu_rhs);
    dxAff = Hbar_diagonal_inverse .* (-rLbar + A*dyAff);

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

    
    rL = g-A*y-(z(1:mIn)-z(mIn+1:mIn*2));
    rA = b-A'*x;
    rC = s+[l; -u] - [x; -x];
    
    % Initial dual gap
    dualGap = (z'*s)/(2*mIn);
    dualGap0 = dualGap;
    
    for i = 1:max_iter
        iter = iter + 1;
        zsl = zl./sl;
        zsu = zu./su;
        
        % Affine step
        rCs = (rC-s);
        rLbar = rL - zsl.*rCs(1:mIn) +zsu.*rCs(1+mIn:2*mIn);
      
        % Calculate the factor in the normal equation       
		Hbar_diagonal_inverse = 1./(zsl+zsu);
		normalfactor = A' * (Hbar_diagonal_inverse .*  A);
        R = chol(normalfactor);
        
        mu_rhs = rA + A' * (Hbar_diagonal_inverse .* rLbar);
		dyAff = R \ (R' \ mu_rhs);
		dxAff = Hbar_diagonal_inverse .* (-rLbar + A*dyAff);
        
        
        dzAff = - [zsl.*dxAff; -zsu.*dxAff] + (z./s).*rCs;
        dsAff = -s-(s./z).*dzAff;
        
        %compute max alpha affine
        dZS = [dzAff; dsAff];
        alphas = (-[z;s]./dZS);
        alphaAff = min([1;alphas(dZS<0)]);
           
        dualGapAff = ((z+alphaAff*dzAff)'*(s+alphaAff*dsAff))/(2*mIn);
        sigma = (dualGapAff/dualGap)^3;
        
        
        % Affine-Centering-Correction Direction
        rSZz = s + dsAff.*dzAff./z-dualGap*sigma*e./z;
        rCs = (rC-rSZz);
        
		rLbar = rL - zsl.*rCs(1:mIn) +zsu.*rCs(1+mIn:2*mIn);
        
        % Calculate the factor in the normal equation
		normalfactor = A' * (Hbar_diagonal_inverse .*  A);
        R = chol(normalfactor);
        
        mu_rhs = rA + A' * (Hbar_diagonal_inverse .* rLbar);
        
        %This is normal equation stuff as well
		dy = R \ (R' \ mu_rhs);
		dx = Hbar_diagonal_inverse .* (-rLbar + A*dy);
        
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


        rL = g-A*y-(z(1:mIn)-z(mIn+1:mIn*2));
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

