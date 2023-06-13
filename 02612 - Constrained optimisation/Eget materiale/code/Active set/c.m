function [val] = c(A,b,x,Wk)
val = A(:,Wk)'*x-b(Wk);
end

