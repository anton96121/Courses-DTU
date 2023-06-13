function [KKT] = get_KKT(H,g,A,b)
% get_KKT   Generate a KKT matrix
%
%
% Syntax: [KKT] = get_KKT(H,g,A,b)
%
%         KKT             : The KKT matrix

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    KKT = [H -A;-A', zeros(size(A,2), size(A,2))];
end