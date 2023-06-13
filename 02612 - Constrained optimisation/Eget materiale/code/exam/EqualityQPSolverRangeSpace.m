function [x, lambda, time_R] = EqualityQPSolverRangeSpace(H,g,A,b)
% EqualityQPSolverRangeSpace   Range Space solver
%
%          min  x'*H*x+g'x
%           x
%          s.t. A x  = b      (Lagrange multiplier: lambda)
%
%
% Syntax: [x, lambda, time_R] = EqualityQPSolverRangeSpace(H,g,A,b)
%
%         x               : Solution
%         lambda          : Lagrange multipier
%         time_R          : Time spend on cholesky factorization

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    % Factorize H
    start = cputime;
    R=chol(H); 
    time_R = cputime-start;
    
    % Solve for x and lambda
    mu=R'\g;
    Hg=R\mu;
    mu=R'\A;
    HA=R\mu;
    lambda = (A'*HA)\(b+A'*Hg);
    x = HA*lambda-Hg;