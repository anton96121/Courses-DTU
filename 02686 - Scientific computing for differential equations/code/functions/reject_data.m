function [hs, idx_reject] = reject_data(hs1)
idx_h = find(hs1.AcceptStep==1);
hs = hs1.h(idx_h);
idx_reject = zeros(1,length(hs1.AcceptStep(hs1.AcceptStep==0)));
i = 1;
last = 1;
l = 1;
for j=1:length(hs1.AcceptStep)
    if last == 1 && hs1.AcceptStep(j) == 0
        idx_reject(l) = i;
        l = l+1;
        last = 0;
    elseif hs1.AcceptStep(j) == 0
        last = 0;
    else
        last = 1;
        i = i+1;
    end
end
idx_reject = idx_reject(1:l-1);
end