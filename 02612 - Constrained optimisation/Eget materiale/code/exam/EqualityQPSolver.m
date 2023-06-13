function [x, lambda,time] = EqualityQPSolver(H, g, A, b, solver)
% EqualityQPSolver      The solver interface for the equality EQP solvers
%
%          min  x'*H*x+g'x
%           x
%          s.t. A x  = b      (Lagrange multiplier: lambda)
%
%
% Syntax: [x, lambda,time] = EqualityQPSolver(H, g, A, b, solver)
%
%         x             : Solution
%         z             : Lagrange multipliers
%         time          : Time used on factorization in some of the
%                           algorithms

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
time = 0;
if solver == "LDLdense"
    [x, lambda] = EqualityQPSolverLDLdense(H,g,A,b);
elseif solver == "LDLsparse"
    [x, lambda] = EqualityQPSolverLDLsparse(H,g,A,b);
elseif solver == "LUdense"
    [x, lambda] = EqualityQPSolverLUdense(H,g,A,b);
elseif solver == "LUsparse"
    [x, lambda] = EqualityQPSolverLUsparse(H,g,A,b);
elseif solver == "rangespace"
    [x, lambda,time] = EqualityQPSolverRangeSpace(H,g,A,b);
elseif solver == "nullspace"
    [x, lambda,time] = EqualityQPSolverNullSpace(H,g,A,b);
else
    error("solver " + solver + "does not exist; possible values are LDLdense, LDLsparse, LUdense, LUsparse, rangespace, and nullspace")
end

