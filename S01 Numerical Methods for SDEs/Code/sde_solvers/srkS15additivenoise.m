function [x,tspan,I01,Z] = srkS15additivenoise(f,g,tspan,x0,R,Z)

  % Number of steps
  steps = numel(tspan);
  
  % Allocate space
  x = zeros(size(x0,1),steps);
  I01 = zeros(size(x0,1),steps);

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
    I1  = (R(:,k)-R(:,k-1));%sqrt(dt)*R(1,k);
    I10 = 1/2*dt*(I1-1/sqrt(3)*sqrt(dt)*Z(1,k)); %
    

    % H values
    H1 = x(:,k-1); % same for g and f
    
    fx = f(H1,tspan(k-1));
    gx = g(H1,tspan(k-1) + 1*dt);

    H02 = x(:,k-1) + fx*dt;
    
    fxx = f(H02,tspan(k-1)+1*dt);
    gxx = g(H02,tspan(k-1));
    H03 = x(:,k-1) + 1/4*(fx+fxx)*dt+(gx+1/2*gxx)*I10/dt;
    
    fxxx = f(H03,tspan(k-1)+1/2*dt);
    
    % Step
    x(:,k) = x(:,k-1) + ...
        (1/6*fx+1/6*fxx+2/3*fxxx)*dt + ...
        (1*I1 + 1*I10)*gx + ... 
        (-1*I10) * gxx;
    I01(:,k) = I10;
    

  end

  

