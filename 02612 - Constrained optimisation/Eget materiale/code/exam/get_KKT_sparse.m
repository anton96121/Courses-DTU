function [KKT] = get_KKT_sparse(H,g,A,b)
% get_KKT_sparse   Generate a sparse KKT matrix
%
%
% Syntax: [KKT] = get_KKT_sparse(H,g,A,b)
%
%         KKT             : The sparse KKT matrix

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    KKT = sparse([H -A;-A', zeros(size(A,2), size(A,2))]);
end