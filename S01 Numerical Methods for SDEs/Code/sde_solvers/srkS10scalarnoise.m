function [x,tspan] = srkS10scalarnoise(f,g,tspan,x0,R)

  % Number of steps
  steps = numel(tspan);
  
  % Allocate space
  x = zeros(size(x0,1),steps);

  % Initial state
  x(:,1) = x0;

  
  % Iterate
  for k=2:steps

    % Time discretization
    dt = tspan(k)-tspan(k-1);

    % Increment
    db  = (R(:,k)-R(:,k-1));
    dbb = 1/2*(db^2 - dt);
        
    % Evaluate only once
    fx = f(x(:,k-1),tspan(k-1));
    gx = g(x(:,k-1),tspan(k-1));
    
    % Supporting values
    x2  = x(:,k-1) + fx*dt;
    tx2 = x2 + gx*dbb/sqrt(dt);
    tx3 = x2 - gx*dbb/sqrt(dt);
    
    % Evaluate the remaining values
    fx2 = f(x2,tspan(k-1)+dt);
    gx2 = g(tx2,tspan(k-1)+dt);
    gx3 = g(tx3,tspan(k-1)+dt);
    
    % Step
    x(:,k) = x(:,k-1) + ...
        (fx+fx2)*dt/2 + ...
        gx*db + ...
        sqrt(dt)/2*(gx2 - gx3);
    
  end

