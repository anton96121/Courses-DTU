function [x,tspan] = srkS15scalarnoise(f,g,tspan,x0,R)

  % Number of steps
  steps = numel(tspan);
  
  % Allocate space
  x = zeros(size(x0,1),steps);

  % Initial state
  x(:,1) = x0;

  if nargin<6
    Z = randn(1,steps);
  end
  
  % Iterate
  for k=2:steps

    % Time discretization
    dt = tspan(k)-tspan(k-1);

    % Increment
    I1  = (R(:,k)-R(:,k-1));
    I11 = 1/2*(I1^2 - dt);
    I111 = 1/6*(I1^3-3*dt*I1);
    I01 = 1/2*dt*(I1-1/sqrt(3)*sqrt(dt)*Z(1,k)); 

    

    % H values
    H1 = x(:,k-1); % same for g and f
    
    fx = f(H1,tspan(k-1));
    gx = g(H1,tspan(k-1));

    H02 = x(:,k-1) + fx*dt;
    H12 = x(:,k-1) + 1/4*fx*dt-1/2*gx*sqrt(dt);
    
    fxx = f(H02,tspan(k-1)+1*dt);
    gxx = g(H12,tspan(k-1)+1/4*dt);
    H03 = x(:,k-1) + 1/4*(fx+fxx)*dt+(gx+1/2*gxx)*I01/dt;
    H13 = x(:,k-1) + fx*dt + gx*sqrt(dt);
    
    fxxx = f(H03,tspan(k-1)+1*dt);
    gxxx = g(H13,tspan(k-1)+dt);
    % H04 = H1
    H14 = x(:,k-1) + 1/4*fxxx*dt + sqrt(dt)*(2*gx-gxx+1/2*gxxx);
    
    % fxxxx = fx
    gxxxx = g(H14,tspan(k-1)+1/4*dt);
    % Step
    x(:,k) = x(:,k-1) + (1/6*fx+1/6*fxx+2/3*fxxx)*dt + (-1*I1+1*I11/sqrt(dt)+2*I01/dt-2*I111/dt)*gx+(4/3*I1-4/3*I11/sqrt(dt)-4/3*I01/dt+5/3*I111/dt)*gxx+(2/3*I1+1/3*I11/sqrt(dt)-2/3*I01/dt-2/3*I111/dt)*gxxx+I111/dt*gxxxx;
    
  end

