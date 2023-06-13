function [x,f_opt,r] = simplex(A,b,c,B,blands,display)
if nargin<5
  blands = true;
end
if nargin<6
  display = true;
end



n = size(A,2);
maxits = 1000;

b = b(:);
c = c(:);
N = setdiff(1:n,B);

for j=1:maxits
    A_B = A(:,B);
    A_N = A(:,N);
    
    [L,U,P]  = lu(A_B);
   
    r = c(N)-A_N'*(L'\(U'\(P*c(B))));
    
    if display
        fprintf('Iteration %2d: B=[',j)
        for i = 1:length(B)-1
            fprintf('%2d,',B(i))
        end
        fprintf('%2d',B(end))
        fprintf('] f_0=%.2f', c(B)'*(U\(L\(P'*b))))
        fprintf(' r=[');
        fprintf('%g, ', r(1:end-1));
        fprintf('%g]\n', r(end));
    end
    
    
    if all(r >= 0)
        if display
            fprintf('Optimal\n')
        end
        break
    end
    
    I = find(r<0);
    [~,idx] = min(r(I));
    idxN = I(idx);
    iEntering = N(idxN);
    
    
    d = -U\(L\(P'*A(:,iEntering)));
    ratios = (U\(L\(P'*b)))./(-d);
    ratios(d>=0) = inf;
    
    if blands && any(ratios == 0)
        iEntering = min(N(r<0));
        d = -U\(L\(P'*A(:,iEntering)));
        ratios = (U\(L\(P'*b)))./(-d);
        ratios(d>=0) = inf;
        
        if any(ratios == 0)
            [iLeaving,idxB] = min(B(ratios==0));
        else
            [~,idxB] = min(ratios);
            iLeaving = B(idxB);
        end
    else
        [~,idxB] = min(ratios);
         iLeaving = B(idxB);
    end
    B(idxB) = iEntering;
    N(idxN) = iLeaving;
end
if ~all(r>=0)
    warning('not at optimum after 1000iters')
end

xB = (U\(L\(P'*b)));
f_opt = c(B)'*xB;
x = zeros(n,1);
x(B) = xB;
end

