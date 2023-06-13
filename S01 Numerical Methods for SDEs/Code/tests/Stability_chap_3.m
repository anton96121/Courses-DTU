
precision = 0.002;
xlimit = [-3,3];
ylimit = [-3,3];
% RK4
A = [0];
b = [1]';

main_title = "Euler-Maruyama";
stability(A,b,xlimit,ylimit,precision, main_title)


%
precision = 0.002;
xlimit = [-3,3];
ylimit = [-3,3];

% SRA3
A = [0 0 0; 1 0 0 ; 1/4 1/4 0];
b = [1/6 1/6 2/3]';

main_title = "SRA3";
stability(A,b,xlimit,ylimit,precision, main_title)

%% bias of the variance
bias_SRA3 = @(mu,sigma,dt) sigma^2*dt^2*mu*(dt^3*mu^3 - 6*dt^2*mu^2 + 21*dt*mu - 24)/(2*dt^5*mu^5 - 12*dt^4*mu^4 + 42*dt^3*mu^3 - 96*dt^2*mu^2 + 144*dt*mu - 144);
bias_EM = @(mu,sigma,dt) dt*sigma^2/(-2*dt*mu - 4);

mu_try = 4;
sigma_try = 100;
dt_vec = [0.0001:0.0001:0.001,0.002:0.001:0.01,0.02:0.01:0.1,0.2:0.1:0.5];
bias_points_EM = zeros(1,length(dt_vec));
bias_points_SRA3 = zeros(1,length(dt_vec));

for i = 1:length(dt_vec)
    bias_points_EM(i) = abs(bias_EM(mu_try,sigma_try,dt_vec(i)));
    bias_points_SRA3(i) = abs(bias_SRA3(mu_try,sigma_try,dt_vec(i))); 
end

figure
loglog(dt_vec,bias_points_SRA3)
hold on 
loglog(dt_vec,bias_points_EM)
xlabel('$\Delta t$','interpreter','latex','FontSize',14); ylabel('Bias','interpreter','latex','FontSize',14)
legend('SRA3','EM','Location','northwest','interpreter','latex','FontSize',12)
title('Bias of the variance', 'Interpreter','latex','FontSize',20)
grid()


p_SRA3 = polyfit(log(dt_vec),log(bias_points_SRA3),1);
p_EM = polyfit(log(dt_vec),log(bias_points_EM),1);
disp(p_EM(1))
disp(p_SRA3(1))

%%
function [] = stability(A,b,xlimit,ylimit, precision, main_title)
    s = length(b);
    re = [xlimit(1):precision:xlimit(2)];
    im = [ylimit(1):precision:ylimit(2)];
    abs_matrix = zeros(length(re),length(im));
    fprintf('\n');
    reverseStr = '';
    for i = 1:length(re)
        for j = 1:length(im)
            abs_matrix(i,j) = abs(1+complex(re(i),im(j))*b'*inv(eye(s)-complex(re(i),im(j))*A)*ones(s,1));

        end
       percentDone = 100 * i / length(re);
       msg = sprintf('Percent done: %3.1f', percentDone);
       fprintf([reverseStr, msg]);
       reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    fprintf('\n');
    plot_matrix = abs_matrix;
    %plot_matrix = (abs_matrix > 1);
    plot_matrix(plot_matrix>1) = 10.0;
    %plot_matrix(round(plot_matrix,1)==1.0) = 0;

    figure
    hold on 
    cmap = winter(256);
    cmap(end,:) = [0.95,0.95,0.95]';
    cmap = [[0 0.95]',[0.4470 0.95]',[0.7410 0.95]'];
    
    colormap(cmap)
    %colorbar;
    stabMap = imagesc(re',im',plot_matrix');
    alpha(stabMap,1);
    grid on
    plot(re,zeros(1,length(re)),'k', 'LineWidth', 0.5)
    plot(zeros(1,length(im)),im,'k', 'LineWidth', 0.5)
    xlim([min(re) max(re)])
    ylim([min(im) max(im)])
    title(main_title, 'Interpreter','latex','FontSize',20)
    xlabel("Re($\lambda$)", 'Interpreter','latex','FontSize',14)
    ylabel("Im($\lambda$)", 'Interpreter','latex','FontSize',14)
    set(gca, 'layer', 'top');


    if all(all(abs_matrix(re<0,:)<1))
        disp('The method is A-stable')
        if all(all(abs_matrix(re>0,:)>1))
            disp('Only Re(z)<0 behaves stable')
        end
        if b'*inv(A)*ones(s,1)-1 < 0.000001
            disp('The method is L-stable')
        end
    else
        disp('The method is not A-stable')
    end

end




