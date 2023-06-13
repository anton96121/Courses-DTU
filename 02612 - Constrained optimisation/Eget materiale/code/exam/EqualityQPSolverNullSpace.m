function [x,lambda,time_N] = EqualityQPSolverNullSpace(H,g,A,b)
% EqualityQPSolverNullSpace   Null Space solver
%
%          min  x'*H*x+g'x
%           x
%          s.t. A x  = b      (Lagrange multiplier: lambda)
%
%
% Syntax: [x,lambda,time_N] = EqualityQPSolverNullSpace(H,g,A,b)
%
%         x               : Solution
%         lambda          : Lagrange multipier
%         time_N          : Time spend on qr factorization

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    [n,m] = size(A);
    
    % Factorize A
    start = cputime;
    [Q,R] = qr(A,'vector');
    time_N = cputime-start;
    
    % Solve for x and lambda
    Qrange = Q(:,1:m);
    Qnull = Q(:,m+1:n);
    R = R(1:m,1:m);
    Y = (R'\b);
    Qnt = Qnull';
    lpre = Qnt*H*Qnull;
    L = chol(lpre);
    mu=L'\(-Qnt*(H*Qrange*Y+g));
    Z=L\mu;
    x = Qrange*Y+Qnull*Z;
    lambda = R\Qrange'*(g+H*x);

    
    