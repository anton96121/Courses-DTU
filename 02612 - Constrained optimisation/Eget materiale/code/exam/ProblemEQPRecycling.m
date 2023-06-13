function [H,g,A,b] = ProblemEQPRecycling(n, uhat, d0)
% ProblemEQPRecycling   Generate an instance of the recycling problem 
%                       described in exercise 1.5
%
%
% Syntax: [H,g,A,b] = ProblemEQPRecycling(n, uhat, d0)
%
%         H               : Hessian
%         g               : Linear term of the objective
%         A               : Matrix of the constraints
%         b               : lhs of the constraints 

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    H = eye(n+1);
    g = uhat*ones(n+1,1);
    A = eye(n);
    sub = eye(n);
    sub = circshift(sub,1,2);
    A = A - sub;
    A = padarray(A,[0 1],0,'post');
    A(end-1,end) = -1;
    A(end-1,end-1) = -1;
    A(end-1,end-2) = 1;
    A = A';
    b = zeros(n,1);
    b(n) = -d0;
end
