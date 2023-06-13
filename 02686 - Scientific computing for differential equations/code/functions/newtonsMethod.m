function [x, xdot, function_calls] = newtonsMethod(f, jac, told, Rterm, h, xguess, tolerance, maxiterations, params, function_calls)
%NEWTONSMETHODODE Does up to maxiterations rounds of newtons method to find
% the next step of the implicit within given tolerance.
% In case of slow convergence, it returns the guess it has reached at that
% point.
% f - the ODE
% jac - the jacobian
% told - the last time step
% Rterm - the term in the residual not encompassed by x_n+1 nor hf(x_n+1)
% dt - The step size
% xguess - The initial starting point (x_n+hf(x_n))
% tolerance - The size of the residual at which we exit
% maxiterations - Maximum allowed number of Newton iterations
% params - Parameters for f and jac
i = 0;
t = told + h;
x = xguess;
xdot = f(t, x, params);
J = jac(t,x,params);
function_calls(1) = function_calls +2;
R = x - xdot*h - Rterm;
I = eye(length(Rterm));
while (i < maxiterations) && (max(abs(R)) > tolerance) %Iteratively improve guess using Newton's method.
    
    %The Jacobian tells us the change of the residual in x, and we then try
    %to improve x using that information.
    dRdx = I - J * h;
    dx = dRdx\R;
    x = x - dx;
    xdot = f(t, x, params);
    %Calculate jacobian and residual again
    J = jac(t,x,params);
    function_calls = function_calls +2;
    R = x - h*xdot - Rterm;
    i = i+ 1;
end
if i==maxiterations && (max(abs(R)) > tolerance)
    disp("Not converging..") %big oh no if we're not converging
end

