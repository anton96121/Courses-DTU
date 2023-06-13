function [H,g,A,b,x,lambda] = generateRandomEQP(n,m)
% generateRandomEQP   Generate a random EQP with n variables and m
%                     constraints
%
%
% Syntax: [H,g,A,b,x,lambda] = generateRandomEQP(n,m)
%
%         x               : Solution
%         lambda          : Lagrange multipier
%         H               : Hessian
%         g               : Linear term of the objective
%         A               : Matrix of the constraints
%         b               : lhs of the constraints 

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    % Create a symmetric and pos def Hessian
    H = rand(n);
    H = 0.5*(H+H') + n*eye(n);

    % Create a full rank A matrix
    A = 10*rand(n); 
    A = 0.5*(A+A')+n*eye(n);
    A = A(:,1:m); 

    % Create a random solution to the system
    x = rand(n,1);
    lambda = rand(m,1);

    % Create the matching g and b
    KKT = [H -A;-A' zeros(m)];
    sol = [x; lambda];
    rhs = KKT*sol;

    g = -rhs(1:n);
    b = -rhs(n+1:n+m);
end







