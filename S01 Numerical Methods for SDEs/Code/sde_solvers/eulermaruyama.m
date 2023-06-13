function [x,tspan] = eulermaruyama(f,g,tspan,x0,R)
  
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
    db = R(:,k)-R(:,k-1);

    % Step
    x(:,k) = x(:,k-1) + ...
        f(x(:,k-1),tspan(k-1))*dt + ...
        g(x(:,k-1),tspan(k-1))*db;
    
  end


