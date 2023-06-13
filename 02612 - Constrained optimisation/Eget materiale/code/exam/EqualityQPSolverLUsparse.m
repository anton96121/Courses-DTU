function [x, lambda] = EqualityQPSolverLUsparse(H,g,A,b)
% EqualityQPSolverLUsparse   Sparse LU solver
%
%          min  x'*H*x+g'x
%           x
%          s.t. A x  = b      (Lagrange multiplier: lambda)
%
%
% Syntax: [x, lambda] = EqualityQPSolverLUsparse(H,g,A,b)
%
%         x               : Solution
%         lambda          : Lagrange multipier

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    % Create a sparse KKT system
    KKT = get_KKT_sparse(H,g,A,b);
    
    % Fatorize the KKT matrix
    [L,U,p] = lu(KKT,'vector');
    
    % Solve for x and lambda
    rhs = [g;b];
    solution(p) = U \ (L \ ( -rhs(p)));
    x = solution(1:size(H,1));
    lambda = solution(size(H,1)+1:size(H,1)+size(b,1));
end
