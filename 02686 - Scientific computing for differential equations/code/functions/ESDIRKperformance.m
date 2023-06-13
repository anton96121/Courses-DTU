function ESDIRKperformance(info, stats)
s = info.nStage;
I = find(stats.r > 1) % failed steps that converged
J = find(stats.Diverged > 0) % diverged steps
K = find(stats.SlowConv > 0) % slow converging steps (reached itermax)
%% Plotting (t_n, h_n) where t_{n+1} = t_n + h_n:
figure
subplot(4,2,[1 3])
semilogy(stats.t, stats.h, '.-'), hold on
semilogy(stats.t(I), stats.h(I), '.r')
semilogy(stats.t(J), stats.h(J), 'xr', 'markersize', 8, 'linewidth', 2)
semilogy(stats.t(K), stats.h(K), 'or', 'markersize', 8, 'linewidth', 2)
%xlabel('t_n')
ylabel('h_n')
%legend('t_{n+1} = t_n + h_n', 'RejStep', 'DivStep', 'SlowCon')
axis tight
%% Plotting (t_n, r_n/tol):
subplot(4,2,[5 7])
semilogy(stats.t, stats.r, '.-'), hold on
semilogy(stats.t(I), stats.r(I), '.r')
plot([min(info.tspan(1)) max(info.tspan(2))], [1 1], '-k') % steps above this line is failed steps that converged
xlabel('t_n')
ylabel('r_n/tol')
%legend('t_{n+1} = t_n + h_n','RejStep')
axis tight
%% Plotting (t_n, iter_n) at each stage:
for i = 1:s-1
    if i > 3, break, end
    J = find(stats.Diverged == i+1); % failed steps that converged
    K = find(stats.SlowConv == i+1); % slow converging steps (reached itermax)
    subplot(4,2,2*i) % iterations in stage i+1
    plot(stats.t, stats.iter(:, i+1), '.-b'), hold on
    plot(stats.t(I), stats.iter(I, i+1), '.r')
    plot(stats.t(J), stats.iter(J, i+1), 'xr', 'markersize', 8, 'linewidth', 2)
    plot(stats.t(K), stats.iter(K, i+1), 'or', 'markersize', 8, 'linewidth', 2)
    if i == s-1, xlabel('t_n'), end
    ylabel(['s_',int2str(i+1)])%,': ',int2str(sum(stats.iter(:, i+1), 1))])
    ylim([0 info.iterMax])
    %axis tight
end
%% Adding text:
g{1} = ['nStep : ',int2str(info.nStep)];
g{2} = ['nFail : ',int2str(info.nFail)];
g{3} = ['nDiv  : ',int2str(info.nDiverge)];
g{4} = ['nSlow : ',int2str(info.nSlowConv)];
h{1} = ['nFun  : ',int2str(info.nFun)];
h{2} = ['nJac  : ',int2str(info.nJac)];
h{3} = ['nLU   : ',int2str(info.nLU)];
h{4} = ['nBack : ',int2str(info.nBack)];
e{1} = ['Method : ',info.Method];
e{2} = ['absTol : ',num2str(info.absTol)];
e{3} = ['relTol : ',num2str(info.relTol)];
for i = 1:s-1
    f{i} = ['nIter_',int2str(i+1),' : ',int2str(sum(stats.iter(:, i+1), 1))];
end
subplot(4,2,8)
%plot(abs(stats.eig(:,1)),'x-k'), hold on, plot(abs(stats.eig(:,2)),'o-k')
axis off
ha(4) = text(-0.1,  0.7, e); ha(3) = text(0.6,  0.7, f);
ha(1) = text(-0.1, -0.4, g); ha(2) = text(0.6, -0.4, h);
set(ha, 'fontname', 'courier');
%maximize % maximizes figure window