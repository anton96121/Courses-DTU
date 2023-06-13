addpath(genpath("..\sde_solvers"))

%
dt = 0.001;
r = 0.1;
sigma = 5;
Tend = 20;
steps = ceil(Tend/dt);
N = 200;

closed_est = zeros(2,N);
opt_exactMLE_est = zeros(2,N);
opt_aproxMLE_est = zeros(2,N);
opt_aproxMLE_estSRA3 = zeros(2,N);

% MLE
a_fe = @(r,sigma) exp(-r * dt);
s_fe = @(r,sigma) sigma/(2*r) * (1 - exp(-2*r*dt));

a_f = @(r,sigma) (1 - 1*r * dt);
s_f = @(r,sigma) sigma*dt;

a_f_add = @(r,sigma) (dt^2*r^2-2*dt*r+2)/2; % MÅ IKKE SLETTES. DET ER DEN DÅRLIGE FRA SDE BOG
s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;

a_f_sra3 = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
s_f_sra3 = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;


for j = 1:N
    % Simulate data
    X = zeros(1,steps);
    T = zeros(1,steps);
    
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end

    % Closed-form maximum likelihood, s240
    r_est = -1/dt * log(sum(X(1:end-1) .* X(2:end)) / ...
        sum(X(1:end-1) .* X(1:end-1)));
    sigma_est = 1/(length(X)-1) * (2*r_est / (1 - exp(-2*r_est*dt))) * ...
        sum((X(2:end) - exp(-r_est*dt) * X(1:end-1)).^2);

    closed_est(1,j) = r_est;%-log(r_est)/(Tend/3);
    closed_est(2,j) = sigma_est;
    
    % Optimized maximum likelihood, s241
    ell = @(p) 0.5 * sum(log(2*pi*s_fe(p(1),p(2))) + 1/s_fe(p(1),p(2)) * (X(2:end) - a_fe(p(1),p(2)) * X(1:end-1)).^2);

    p_e = fminsearch(ell,[2 2]);
    opt_exactMLE_est(1,j) = p_e(1);
    opt_exactMLE_est(2,j) = p_e(2);

    % Optimized maximum likelihood from Euler-Maruyama, s244   
    ell_em = @(p) 0.5 * sum(log(2*pi*s_f(p(1),p(2))) + 1/s_f(p(1),p(2)) * (X(2:end) - a_f(p(1),p(2)) * X(1:end-1)).^2);
    %ell_em = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_em = fminsearch(ell_em,[2 2]);
    opt_aproxMLE_est(1,j) = p_em(1);
    opt_aproxMLE_est(2,j) = p_em(2);

    % Optimized maximum likelihood from SRA3  
    %ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_sra3(p(1),p(2))) + 1/s_f_sra3(p(1),p(2)) * (X(2:end) - a_f_sra3(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA3(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA3(2,j) = p_SRA3(2);
end    
%
%
meanClosedForm = mean(closed_est,2);
meanExactMLE = mean(opt_exactMLE_est,2);
meanEM = mean(opt_aproxMLE_est,2);
meanSR3 = mean(opt_aproxMLE_estSRA3,2);

[f,xi]=ksdensity(opt_aproxMLE_estSRA3(1,:));
[fmax, idx] = max(f);
r_SRA3 = xi(idx);

max_mu = max([max(closed_est(1,:)),max(opt_aproxMLE_estSRA3(1,:)),max(opt_aproxMLE_est(1,:)),max(opt_exactMLE_est(1,:))]);
min_mu = min([min(closed_est(1,:)),min(opt_aproxMLE_estSRA3(1,:)),min(opt_aproxMLE_est(1,:)),min(opt_exactMLE_est(1,:))]);

max_sigma = max([max(closed_est(2,:)),max(opt_aproxMLE_estSRA3(2,:)),max(opt_aproxMLE_est(2,:)),max(opt_exactMLE_est(2,:))]);
min_sigma = min([min(closed_est(2,:)),min(opt_aproxMLE_estSRA3(2,:)),min(opt_aproxMLE_est(2,:)),min(opt_exactMLE_est(2,:))]);

mu_delta = max_mu-min_mu;
max_mu = max_mu + mu_delta*0.2;
min_mu = min_mu - mu_delta*0.2;

sigma_delta = max_sigma-min_sigma;
max_sigma = max_sigma + sigma_delta*0.2;
min_sigma = min_sigma - sigma_delta*0.2;

%
figure;
subplot(2, 2, 1); hold on;
plot(closed_est(1,:),closed_est(2,:),'.')
plot(meanClosedForm(1),meanClosedForm(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
xlim([min_mu,max_mu])
ylim([min_sigma, max_sigma])
title("Closed form",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

subplot(2, 2, 2); hold on;
plot(opt_exactMLE_est(1,:),opt_exactMLE_est(2,:),'.')
plot(meanExactMLE(1),meanExactMLE(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
xlim([min_mu,max_mu])
ylim([min_sigma, max_sigma])
title("Exact MLE",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

subplot(2, 2, 3); hold on;
plot(opt_aproxMLE_est(1,:),opt_aproxMLE_est(2,:),'.')
plot(meanEM(1),meanEM(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
xlim([min_mu,max_mu])
ylim([min_sigma, max_sigma])
title("EM",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

subplot(2, 2, 4); hold on;
plot(opt_aproxMLE_estSRA3(1,:),opt_aproxMLE_estSRA3(2,:),'.')
plot(meanSR3(1),meanSR3(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
xlim([min_mu,max_mu])
ylim([min_sigma, max_sigma])
title("SRA3",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

% add a bit space to the figure
fig = gcf;
% add legend
Lgnd = legend("","Mean estimate", "True parameters",'interpreter','latex');
Lgnd.FontSize = 12;
hold off;


%


figure; subplot(2,3,1);hold on;
plot(closed_est(1,:),closed_est(2,:),'.')
plot(meanClosedForm(1),meanClosedForm(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
xlim([min_mu,max_mu])
ylim([min_sigma, max_sigma])
title("Closed form",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')
subplot(2,3,2);hold on;
plot(opt_exactMLE_est(1,:),opt_exactMLE_est(2,:),'.')
plot(meanExactMLE(1),meanExactMLE(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
xlim([min_mu,max_mu])
ylim([min_sigma, max_sigma])
title("Exact MLE",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma$",'interpreter','latex')
%legend("","Mean estimate", "True parameters",'Location','southeast','interpreter','latex')
subplot(2,3,4);hold on;
plot(opt_aproxMLE_est(1,:),opt_aproxMLE_est(2,:),'.')
plot(meanEM(1),meanEM(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
xlim([min_mu,max_mu])
ylim([min_sigma, max_sigma])
title("EM",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma$",'interpreter','latex')
%legend("","Mean estimate", "True parameters",'Location','southeast','interpreter','latex')
subplot(2,3,5);hold on;
plot(opt_aproxMLE_estSRA3(1,:),opt_aproxMLE_estSRA3(2,:),'.')
plot(meanSR3(1),meanSR3(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
xlim([min_mu,max_mu])
ylim([min_sigma, max_sigma])
title("SRA3",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma$",'interpreter','latex')
%legend("","Mean estimate", "True parameters",'Location','southeast','interpreter','latex')

% Create a tile on the right column to get its position
ax = subplot(2,3,3,'Visible','off');
axPos = ax.Position;
delete(ax)
% Construct a Legend with the data from the sub-plots
hL = legend("","Mean estimate", "True parameters",'interpreter','latex');
% Move the legend to the position of the extra axes
hL.Position(1:2) = axPos(1:2);


hold off;

%


figure; hold on;
plot(meanClosedForm(1),meanClosedForm(2),'.','MarkerSize',20)
plot(meanExactMLE(1),meanExactMLE(2),'.','MarkerSize',20)
plot(meanEM(1),meanEM(2),'.','MarkerSize',20)
plot(meanSR3(1),meanSR3(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
max_r = max([abs(meanClosedForm(1)-r),abs(meanExactMLE(1)-r),abs(meanEM(1)-r),abs(meanSR3(1)-r)]);
max_sigma = max([abs(meanClosedForm(2)-sigma),abs(meanExactMLE(2)-sigma),abs(meanEM(2)-sigma),abs(meanSR3(2)-sigma)]);
maxmax = max(max_r,max_sigma);
xlim([r-2*maxmax,r+2*maxmax])
ylim([sigma-2*maxmax,sigma+2*maxmax])
legend("Closed form", "Exact MLE", "Euler-Maruyama", "SRA3", "True parameters",'Location','southeast','interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma$",'interpreter','latex')

%%
Tend_vec = 2:2:8;
SRA3_params = zeros(2,length(Tend_vec));

for i = 1:length(Tend_vec)

dt = 0.0001;
r = 2;
sigma = 4;
Tend = Tend_vec(i);
steps = ceil(Tend/dt);
N = 200;

opt_aproxMLE_estSRA3 = zeros(2,N);

a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;

for j = 1:N
    % Simulate data
    X = zeros(1,steps);
    T = zeros(1,steps);
    
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end

    % Optimized maximum likelihood from SRA3  
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA3(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA3(2,j) = p_SRA3(2);
end    
%
%
SRA3_params(:,i) = mean(opt_aproxMLE_estSRA3,2);
disp(i)
end

figure
plot(Tend_vec, (SRA3_params(1,:)-2))
xlim([-5,400])

%% dt fejl




dt_vec = [0.0001:0.0001:0.001,0.002:0.001:0.01,0.02:0.01:0.1,0.2:0.1:1];
SRA3_params = zeros(2,length(dt_vec));
Closed_params = zeros(2,length(dt_vec));
MLE_params = zeros(2,length(dt_vec));
EM_params = zeros(2,length(dt_vec));

for i = 1:length(dt_vec)

dt = dt_vec(i);
r = 2;
sigma = 2;
Tend = 200;
steps = ceil(Tend/dt);
N = 200;

closed_est = zeros(2,N);
opt_exactMLE_est = zeros(2,N);
opt_aproxMLE_est = zeros(2,N);
opt_aproxMLE_estSRA3 = zeros(2,N);

% MLE
a_fe = @(r,sigma) exp(-r * dt);
s_fe = @(r,sigma) sigma/(2*r) * (1 - exp(-2*r*dt));

a_f = @(r,sigma) (1 - 1*r * dt);
s_f = @(r,sigma) sigma*dt;

a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;


for j = 1:N
    % Simulate data
    X = zeros(1,steps);
    T = zeros(1,steps);
    
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end

    % Closed-form maximum likelihood, s240
    r_est = -1/dt * log(sum(X(1:end-1) .* X(2:end)) / ...
        sum(X(1:end-1) .* X(1:end-1)));
    sigma_est = 1/(length(X)-1) * (2*r_est / (1 - exp(-2*r_est*dt))) * ...
        sum((X(2:end) - exp(-r_est*dt) * X(1:end-1)).^2);

    closed_est(1,j) = r_est;
    closed_est(2,j) = sigma_est;
    
    % Optimized maximum likelihood, s241
    ell = @(p) 0.5 * sum(log(2*pi*s_fe(p(1),p(2))) + 1/s_fe(p(1),p(2)) * (X(2:end) - a_fe(p(1),p(2)) * X(1:end-1)).^2);

    p_e = fminsearch(ell,[2 2]);
    opt_exactMLE_est(1,j) = p_e(1);
    opt_exactMLE_est(2,j) = p_e(2);

    % Optimized maximum likelihood from Euler-Maruyama, s244   
    ell_em = @(p) 0.5 * sum(log(2*pi*s_f(p(1),p(2))) + 1/s_f(p(1),p(2)) * (X(2:end) - a_f(p(1),p(2)) * X(1:end-1)).^2);

    p_em = fminsearch(ell_em,[2 2]);
    opt_aproxMLE_est(1,j) = p_em(1);
    opt_aproxMLE_est(2,j) = p_em(2);

    % Optimized maximum likelihood from SRA3  
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA3(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA3(2,j) = p_SRA3(2);

end    
%
%
Closed_params(:,i) = mean(closed_est,2);
MLE_params(:,i) = mean(opt_exactMLE_est,2);
EM_params(:,i) = mean(opt_aproxMLE_est,2);
SRA3_params(:,i) = mean(opt_aproxMLE_estSRA3,2);
%
%
disp(i)
end

%%
figure
semilogx(dt_vec, Closed_params(1,:), "Color","#0072BD", 'LineWidth',1.5)
hold on
semilogx(dt_vec, Closed_params(1,:),'.', "Color","#0072BD", ...
    'MarkerSize',20)
semilogx(dt_vec, MLE_params(1,:), "Color","#D95319", 'LineWidth',1.5)
semilogx(dt_vec, MLE_params(1,:),'.', "Color","#D95319", ...
    'MarkerSize',20)
semilogx(dt_vec, EM_params(1,:), "Color","#7E2F8E", 'LineWidth',1.5)
semilogx(dt_vec, EM_params(1,:),'.', "Color","#7E2F8E", ...
    'MarkerSize',20)
semilogx(dt_vec, SRA3_params(1,:), "Color","#77AC30", 'LineWidth',1.5)
semilogx(dt_vec, SRA3_params(1,:),'.', "Color","#77AC30", ...
    'MarkerSize',20)
yline(2,'--',"Color","#4DBEEE", 'LineWidth',2)
ylim([1.7,2.3])
grid()
ylabel("$\mu$",'interpreter','latex')
xlabel("$\Delta t$",'interpreter','latex')
legend("","Closed Form","","Exact MLE","","EM MLE", "","SRA3 MLE","True parameter",'interpreter','latex');


figure
semilogx(dt_vec, Closed_params(1,:), "Color","#0072BD", 'LineWidth',1.5)
hold on
semilogx(dt_vec, MLE_params(1,:), "Color","#D95319", 'LineWidth',1.5)
semilogx(dt_vec, EM_params(1,:), "Color","#7E2F8E", 'LineWidth',1.5)
semilogx(dt_vec, SRA3_params(1,:), "Color","#77AC30", 'LineWidth',1.5)
yline(2,'--',"Color","#4DBEEE", 'LineWidth',2)
ylim([1.8,2.2])
xlim([0.0001,1.3])
grid()
ylabel("$\mu$",'interpreter','latex')
xlabel("$\Delta t$",'interpreter','latex')
legend("Closed Form","Exact MLE","EM MLE","SRA3 MLE","True parameter",'interpreter','latex', 'location', 'northwest');

figure
semilogx(dt_vec, Closed_params(2,:), "Color","#0072BD", 'LineWidth',1.5)
hold on
semilogx(dt_vec, MLE_params(2,:), "Color","#D95319", 'LineWidth',1.5)
semilogx(dt_vec, EM_params(2,:), "Color","#7E2F8E", 'LineWidth',1.5)
semilogx(dt_vec, SRA3_params(2,:), "Color","#77AC30", 'LineWidth',1.5)
yline(2,'--',"Color","#4DBEEE", 'LineWidth',2)
ylim([1.8,2.2])
xlim([0.0001,1.3])
grid()
ylabel("$\sigma^2$",'interpreter','latex')
xlabel("$\Delta t$",'interpreter','latex')
legend("Closed Form","Exact MLE","EM MLE","SRA3 MLE","True parameter",'interpreter','latex', 'location', 'northwest');

%%
dt = 0.01;
r = 2;
sigma = 4;
Tend = 500;
steps = ceil(Tend/dt);
N = 200;

opt_aproxMLE_estSRA310 = zeros(2,N);
opt_aproxMLE_estSRA350 = zeros(2,N);
opt_aproxMLE_estSRA3250 = zeros(2,N);
opt_aproxMLE_estSRA31000 = zeros(2,N);

a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;

for j = 1:N
    % Simulate data
    X = zeros(1,steps);
    T = zeros(1,steps);
    
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end

    % Optimized maximum likelihood from SRA3  
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA31000(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA31000(2,j) = p_SRA3(2);

    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:250/dt) - a_f_add(p(1),p(2)) * X(1:250/dt-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA3250(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA3250(2,j) = p_SRA3(2);

    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:50/dt) - a_f_add(p(1),p(2)) * X(1:50/dt-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA350(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA350(2,j) = p_SRA3(2);

    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:10/dt) - a_f_add(p(1),p(2)) * X(1:10/dt-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA310(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA310(2,j) = p_SRA3(2);
end    


meanSRA310 = mean(opt_aproxMLE_estSRA310,2);
meanSRA350 = mean(opt_aproxMLE_estSRA350,2);
meanSRA3250 = mean(opt_aproxMLE_estSRA3250,2);
meanSRA31000 = mean(opt_aproxMLE_estSRA31000,2);

[f,xi]=ksdensity(opt_aproxMLE_estSRA3(1,:));
[fmax, idx] = max(f);
r_SRA3 = xi(idx);

max_mu = max([max(opt_aproxMLE_estSRA310(1,:)),max(opt_aproxMLE_estSRA350(1,:)),max(opt_aproxMLE_estSRA3250(1,:)),max(opt_aproxMLE_estSRA31000(1,:))]);
min_mu = min([min(opt_aproxMLE_estSRA310(1,:)),min(opt_aproxMLE_estSRA350(1,:)),min(opt_aproxMLE_estSRA3250(1,:)),min(opt_aproxMLE_estSRA31000(1,:))]);

max_sigma = max([max(opt_aproxMLE_estSRA310(2,:)),max(opt_aproxMLE_estSRA350(2,:)),max(opt_aproxMLE_estSRA3250(2,:)),max(opt_aproxMLE_estSRA31000(2,:))]);
min_sigma = min([min(opt_aproxMLE_estSRA310(2,:)),min(opt_aproxMLE_estSRA350(2,:)),min(opt_aproxMLE_estSRA3250(2,:)),min(opt_aproxMLE_estSRA31000(2,:))]);

mu_delta = max_mu-min_mu;
max_mu = max_mu + mu_delta*0.2;
min_mu = min_mu - mu_delta*0.2;

sigma_delta = max_sigma-min_sigma;
max_sigma = max_sigma + sigma_delta*0.2;
min_sigma = min_sigma - sigma_delta*0.2;

%%
figure;
subplot(2, 2, 1); hold on;
plot(opt_aproxMLE_estSRA310(1,:),opt_aproxMLE_estSRA310(2,:),'.')
plot(meanSRA310(1),meanSRA310(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
grid()
title("T = 10",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

subplot(2, 2, 2); hold on;
plot(opt_aproxMLE_estSRA350(1,:),opt_aproxMLE_estSRA350(2,:),'.')
plot(meanSRA350(1),meanSRA350(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
grid()
title("T = 50",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

subplot(2, 2, 3); hold on;
plot(opt_aproxMLE_estSRA3250(1,:),opt_aproxMLE_estSRA3250(2,:),'.')
plot(meanSRA3250(1),meanSRA3250(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
grid()
title("T = 250",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

% add a bit space to the figure
fig = gcf;
% add legend
Lgnd = legend("","", "True parameters",'interpreter','latex','Orientation','horizontal','Box','off');
Lgnd.FontSize = 10;

subplot(2, 2, 4); hold on;
plot(opt_aproxMLE_estSRA31000(1,:),opt_aproxMLE_estSRA31000(2,:),'.')
plot(meanSRA31000(1),meanSRA31000(2),'.','MarkerSize',20)
plot(r,sigma,'.','MarkerSize',20)
grid()
title("T = 500",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

% add a bit space to the figure
fig = gcf;
% add legend
Lgnd = legend("","Mean estimate", "",'interpreter','latex','Orientation','horizontal','Box','off');
Lgnd.FontSize = 10;
hold off;
%

figure;
hold on;
p1 = scatter(opt_aproxMLE_estSRA310(1,:),opt_aproxMLE_estSRA310(2,:),5,'filled', "MarkerEdgeColor","#0072BD", ...
    "MarkerFaceColor","#0072BD");
plot(meanSRA310(1),meanSRA310(2),'.', "MarkerEdgeColor","#0072BD", ...
    "MarkerFaceColor","#0072BD",'MarkerSize',20)
p2 = scatter(opt_aproxMLE_estSRA350(1,:),opt_aproxMLE_estSRA350(2,:),5,'filled', "MarkerEdgeColor","#D95319", ...
    "MarkerFaceColor","#D95319");
plot(meanSRA350(1),meanSRA350(2),'.', "MarkerEdgeColor","#D95319", ...
    "MarkerFaceColor","#D95319",'MarkerSize',20)
p3 = scatter(opt_aproxMLE_estSRA3250(1,:),opt_aproxMLE_estSRA3250(2,:),5,'filled', "MarkerEdgeColor","#7E2F8E", ...
    "MarkerFaceColor","#7E2F8E");
plot(meanSRA3250(1),meanSRA3250(2),'.', "MarkerEdgeColor","#7E2F8E", ...
    "MarkerFaceColor","#7E2F8E",'MarkerSize',20)
p4 = scatter(opt_aproxMLE_estSRA31000(1,:),opt_aproxMLE_estSRA31000(2,:),5,'filled', "MarkerEdgeColor","#77AC30", ...
    "MarkerFaceColor","#77AC30");
plot(meanSRA31000(1),meanSRA31000(2),'.', "MarkerEdgeColor","#77AC30", ...
    "MarkerFaceColor","#77AC30",'MarkerSize',20)
set(p1,'MarkerFaceAlpha',.2);
set(p1,'MarkerEdgeAlpha',.2);
set(p2,'MarkerFaceAlpha',.2);
set(p2,'MarkerEdgeAlpha',.2);
set(p3,'MarkerFaceAlpha',.2);
set(p3,'MarkerEdgeAlpha',.2);
set(p4,'MarkerFaceAlpha',.2);
set(p4,'MarkerEdgeAlpha',.2);
plot(r,sigma,'.','MarkerSize',20, "MarkerEdgeColor","#4DBEEE", ...
    "MarkerFaceColor","#4DBEEE")
grid()
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')
ylim([3.8,4.2])
xlim([1.95,2.15])
Lgnd = legend("","T = 10","", "T = 50","","T = 250","","T = 500","True",'interpreter','latex');

%%

%%
dt = 0.01;
r = 2;
sigma = 4;
Tend = 20;
steps = ceil(Tend/dt);
N = 200;

opt_aproxMLE_estSRA310 = zeros(2,N);
opt_aproxMLE_estSRA350 = zeros(2,N);
opt_aproxMLE_estSRA3250 = zeros(2,N);
opt_aproxMLE_estSRA31000 = zeros(2,N);

a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;

for j = 1:N
    % Simulate data
    X = zeros(1,steps);
    T = zeros(1,steps);
    sigma = 0.1;
    a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
    s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end

    % Optimized maximum likelihood from SRA3  
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA31000(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA31000(2,j) = p_SRA3(2);


    X = zeros(1,steps);
    T = zeros(1,steps);
    sigma = 1;
    a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
    s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA3250(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA3250(2,j) = p_SRA3(2);

    X = zeros(1,steps);
    T = zeros(1,steps);
    sigma = 2.5;
    a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
    s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end    
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA350(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA350(2,j) = p_SRA3(2);

    X = zeros(1,steps);
    T = zeros(1,steps);
    sigma = 5;
    a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
    s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end    
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA310(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA310(2,j) = p_SRA3(2);
end    


meanSRA310 = mean(opt_aproxMLE_estSRA310,2);
meanSRA350 = mean(opt_aproxMLE_estSRA350,2);
meanSRA3250 = mean(opt_aproxMLE_estSRA3250,2);
meanSRA31000 = mean(opt_aproxMLE_estSRA31000,2);

[f,xi]=ksdensity(opt_aproxMLE_estSRA3(1,:));
[fmax, idx] = max(f);
r_SRA3 = xi(idx);

max_mu = max([max(opt_aproxMLE_estSRA310(1,:)),max(opt_aproxMLE_estSRA350(1,:)),max(opt_aproxMLE_estSRA3250(1,:)),max(opt_aproxMLE_estSRA31000(1,:))]);
min_mu = min([min(opt_aproxMLE_estSRA310(1,:)),min(opt_aproxMLE_estSRA350(1,:)),min(opt_aproxMLE_estSRA3250(1,:)),min(opt_aproxMLE_estSRA31000(1,:))]);

max_sigma = max([max(opt_aproxMLE_estSRA310(2,:)),max(opt_aproxMLE_estSRA350(2,:)),max(opt_aproxMLE_estSRA3250(2,:)),max(opt_aproxMLE_estSRA31000(2,:))]);
min_sigma = min([min(opt_aproxMLE_estSRA310(2,:)),min(opt_aproxMLE_estSRA350(2,:)),min(opt_aproxMLE_estSRA3250(2,:)),min(opt_aproxMLE_estSRA31000(2,:))]);

mu_delta = max_mu-min_mu;
max_mu = max_mu + mu_delta*0.2;
min_mu = min_mu - mu_delta*0.2;

sigma_delta = max_sigma-min_sigma;
max_sigma = max_sigma + sigma_delta*0.2;
min_sigma = min_sigma - sigma_delta*0.2;

%%
figure;
subplot(2, 2, 1); hold on;
plot(opt_aproxMLE_estSRA31000(1,:),opt_aproxMLE_estSRA31000(2,:),'.')
plot(meanSRA31000(1),meanSRA31000(2),'.','MarkerSize',20)
plot(r,0.1,'.','MarkerSize',20)
title("$\sigma^2$ = 0.1",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

subplot(2, 2, 2); hold on;
plot(opt_aproxMLE_estSRA3250(1,:),opt_aproxMLE_estSRA3250(2,:),'.')
plot(meanSRA3250(1),meanSRA3250(2),'.','MarkerSize',20)
plot(r,1.0,'.','MarkerSize',20)
title("$\sigma^2$ = 1.0",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

subplot(2, 2, 3); hold on;
plot(opt_aproxMLE_estSRA350(1,:),opt_aproxMLE_estSRA350(2,:),'.')
plot(meanSRA350(1),meanSRA350(2),'.','MarkerSize',20)
plot(r,2.5,'.','MarkerSize',20)
title("$\sigma^2$ = 2.5",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

% add a bit space to the figure
fig = gcf;
% add legend
Lgnd = legend("","", "True parameters",'interpreter','latex','Orientation','horizontal','Box','off');
Lgnd.FontSize = 10;

subplot(2, 2, 4); hold on;
plot(opt_aproxMLE_estSRA310(1,:),opt_aproxMLE_estSRA310(2,:),'.')
plot(meanSRA310(1),meanSRA310(2),'.','MarkerSize',20)
plot(r,5.0,'.','MarkerSize',20)
title("$\sigma^2$ = 5.0",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

% add a bit space to the figure
fig = gcf;
% add legend
Lgnd = legend("","Mean estimate", "",'interpreter','latex','Orientation','horizontal','Box','off');
Lgnd.FontSize = 10;
hold off;
%

figure;
hold on;
p1 = scatter(opt_aproxMLE_estSRA31000(1,:),opt_aproxMLE_estSRA31000(2,:),5,'filled', "MarkerEdgeColor","#0072BD", ...
    "MarkerFaceColor","#0072BD");
plot(meanSRA31000(1),meanSRA31000(2),'.', "MarkerEdgeColor","#0072BD", ...
    "MarkerFaceColor","#0072BD",'MarkerSize',20)
xline(meanSRA31000(1), 'Color', "#0072BD")
p2 = scatter(opt_aproxMLE_estSRA3250(1,:),opt_aproxMLE_estSRA3250(2,:),5,'filled', "MarkerEdgeColor","#D95319", ...
    "MarkerFaceColor","#D95319");
plot(meanSRA3250(1),meanSRA3250(2),'.', "MarkerEdgeColor","#D95319", ...
    "MarkerFaceColor","#D95319",'MarkerSize',20)
xline(meanSRA3250(1), 'Color', "#D95319")
p3 = scatter(opt_aproxMLE_estSRA350(1,:),opt_aproxMLE_estSRA350(2,:),5,'filled', "MarkerEdgeColor","#7E2F8E", ...
    "MarkerFaceColor","#7E2F8E");
plot(meanSRA350(1),meanSRA350(2),'.', "MarkerEdgeColor","#7E2F8E", ...
    "MarkerFaceColor","#7E2F8E",'MarkerSize',20)
xline(meanSRA350(1), 'Color', "#7E2F8E")
p4 = scatter(opt_aproxMLE_estSRA310(1,:),opt_aproxMLE_estSRA310(2,:),5,'filled', "MarkerEdgeColor","#77AC30", ...
    "MarkerFaceColor","#77AC30");
plot(meanSRA310(1),meanSRA310(2),'.', "MarkerEdgeColor","#77AC30", ...
    "MarkerFaceColor","#77AC30",'MarkerSize',20)
xline(meanSRA310(1), 'Color', "#77AC30")
set(p1,'MarkerFaceAlpha',.2);
set(p1,'MarkerEdgeAlpha',.2);
set(p2,'MarkerFaceAlpha',.2);
set(p2,'MarkerEdgeAlpha',.2);
set(p3,'MarkerFaceAlpha',.2);
set(p3,'MarkerEdgeAlpha',.2);
set(p4,'MarkerFaceAlpha',.2);
set(p4,'MarkerEdgeAlpha',.2);
xline(r,'--', 'Color', "#4DBEEE",'LineWidth', 2)
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')
xlim([1.95,2.2])
Lgnd = legend("","$\sigma^2$ = 0.1","","", "$\sigma^2$ = 1.0","","","$\sigma^2$ = 2.5","","","$\sigma^2$ = 5.0","","True",'interpreter','latex');

%%
dt = 0.01;
r = 2;
sigma = 4;
Tend = 20;
steps = ceil(Tend/dt);
N = 200;

opt_aproxMLE_estSRA310 = zeros(2,N);
opt_aproxMLE_estSRA350 = zeros(2,N);
opt_aproxMLE_estSRA3250 = zeros(2,N);
opt_aproxMLE_estSRA31000 = zeros(2,N);

a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;

for j = 1:N
    % Simulate data
    X = zeros(1,steps);
    T = zeros(1,steps);
    r = 0.1;
    a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
    s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end

    % Optimized maximum likelihood from SRA3  
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA31000(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA31000(2,j) = p_SRA3(2);


    X = zeros(1,steps);
    T = zeros(1,steps);
    r = 1;
    a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
    s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA3250(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA3250(2,j) = p_SRA3(2);

    X = zeros(1,steps);
    T = zeros(1,steps);
    r = 4;
    a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
    s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end    
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA350(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA350(2,j) = p_SRA3(2);

    X = zeros(1,steps);
    T = zeros(1,steps);
    r = 10;
    a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
    s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;
    x = 5;
    t = 0;
    r_sim = r;
    A = exp(-r_sim*dt);
    Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Sigma) * randn;
        t = t + dt;
    end    
    ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (X(2:end) - a_f_add(p(1),p(2)) * X(1:end-1)).^2);

    p_SRA3 = fminsearch(ell_SRA3,[2 2]);
    opt_aproxMLE_estSRA310(1,j) = p_SRA3(1);
    opt_aproxMLE_estSRA310(2,j) = p_SRA3(2);
end    


meanSRA310 = mean(opt_aproxMLE_estSRA310,2);
meanSRA350 = mean(opt_aproxMLE_estSRA350,2);
meanSRA3250 = mean(opt_aproxMLE_estSRA3250,2);
meanSRA31000 = mean(opt_aproxMLE_estSRA31000,2);

[f,xi]=ksdensity(opt_aproxMLE_estSRA3(1,:));
[fmax, idx] = max(f);
r_SRA3 = xi(idx);

max_mu = max([max(opt_aproxMLE_estSRA310(1,:)),max(opt_aproxMLE_estSRA350(1,:)),max(opt_aproxMLE_estSRA3250(1,:)),max(opt_aproxMLE_estSRA31000(1,:))]);
min_mu = min([min(opt_aproxMLE_estSRA310(1,:)),min(opt_aproxMLE_estSRA350(1,:)),min(opt_aproxMLE_estSRA3250(1,:)),min(opt_aproxMLE_estSRA31000(1,:))]);

max_sigma = max([max(opt_aproxMLE_estSRA310(2,:)),max(opt_aproxMLE_estSRA350(2,:)),max(opt_aproxMLE_estSRA3250(2,:)),max(opt_aproxMLE_estSRA31000(2,:))]);
min_sigma = min([min(opt_aproxMLE_estSRA310(2,:)),min(opt_aproxMLE_estSRA350(2,:)),min(opt_aproxMLE_estSRA3250(2,:)),min(opt_aproxMLE_estSRA31000(2,:))]);

mu_delta = max_mu-min_mu;
max_mu = max_mu + mu_delta*0.2;
min_mu = min_mu - mu_delta*0.2;

sigma_delta = max_sigma-min_sigma;
max_sigma = max_sigma + sigma_delta*0.2;
min_sigma = min_sigma - sigma_delta*0.2;

%%
figure;
subplot(2, 2, 1); hold on;
plot(opt_aproxMLE_estSRA31000(1,:),opt_aproxMLE_estSRA31000(2,:),'.')
plot(meanSRA31000(1),meanSRA31000(2),'.','MarkerSize',20)
plot(0.1,sigma,'.','MarkerSize',20)
title("$\mu$ = 0.1",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

subplot(2, 2, 2); hold on;
plot(opt_aproxMLE_estSRA3250(1,:),opt_aproxMLE_estSRA3250(2,:),'.')
plot(meanSRA3250(1),meanSRA3250(2),'.','MarkerSize',20)
plot(1,sigma,'.','MarkerSize',20)
title("$\mu$ = 1.0",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

subplot(2, 2, 3); hold on;
plot(opt_aproxMLE_estSRA350(1,:),opt_aproxMLE_estSRA350(2,:),'.')
plot(meanSRA350(1),meanSRA350(2),'.','MarkerSize',20)
plot(4.0,sigma,'.','MarkerSize',20)
title("$\mu$ = 4.0",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

% add a bit space to the figure
fig = gcf;
% add legend
Lgnd = legend("","", "True parameters",'interpreter','latex','Orientation','horizontal','Box','off');
Lgnd.FontSize = 10;

subplot(2, 2, 4); hold on;
plot(opt_aproxMLE_estSRA310(1,:),opt_aproxMLE_estSRA310(2,:),'.')
plot(meanSRA310(1),meanSRA310(2),'.','MarkerSize',20)
plot(10.0,sigma,'.','MarkerSize',20)
title("$\mu$ = 10.0",'interpreter','latex')
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')

% add a bit space to the figure
fig = gcf;
% add legend
Lgnd = legend("","Mean estimate", "",'interpreter','latex','Orientation','horizontal','Box','off');
Lgnd.FontSize = 10;
hold off;
%

figure;
hold on;
p1 = scatter(opt_aproxMLE_estSRA31000(1,:)/0.1,opt_aproxMLE_estSRA31000(2,:),5,'filled', "MarkerEdgeColor","#0072BD", ...
    "MarkerFaceColor","#0072BD");
plot(meanSRA31000(1)/0.1,meanSRA31000(2),'.', "MarkerEdgeColor","#0072BD", ...
    "MarkerFaceColor","#0072BD",'MarkerSize',20)
p2 = scatter(opt_aproxMLE_estSRA3250(1,:)/1,opt_aproxMLE_estSRA3250(2,:),5,'filled', "MarkerEdgeColor","#D95319", ...
    "MarkerFaceColor","#D95319");
plot(meanSRA3250(1)/1,meanSRA3250(2),'.', "MarkerEdgeColor","#D95319", ...
    "MarkerFaceColor","#D95319",'MarkerSize',20)
p3 = scatter(opt_aproxMLE_estSRA350(1,:)/4,opt_aproxMLE_estSRA350(2,:),5,'filled', "MarkerEdgeColor","#7E2F8E", ...
    "MarkerFaceColor","#7E2F8E");
plot(meanSRA350(1)/4,meanSRA350(2),'.', "MarkerEdgeColor","#7E2F8E", ...
    "MarkerFaceColor","#7E2F8E",'MarkerSize',20)
p4 = scatter(opt_aproxMLE_estSRA310(1,:)/10,opt_aproxMLE_estSRA310(2,:),5,'filled', "MarkerEdgeColor","#77AC30", ...
    "MarkerFaceColor","#77AC30");
plot(meanSRA310(1)/10,meanSRA310(2),'.', "MarkerEdgeColor","#77AC30", ...
    "MarkerFaceColor","#77AC30",'MarkerSize',20)
set(p1,'MarkerFaceAlpha',.2);
set(p1,'MarkerEdgeAlpha',.2);
set(p2,'MarkerFaceAlpha',.2);
set(p2,'MarkerEdgeAlpha',.2);
set(p3,'MarkerFaceAlpha',.2);
set(p3,'MarkerEdgeAlpha',.2);
set(p4,'MarkerFaceAlpha',.2);
set(p4,'MarkerEdgeAlpha',.2);
plot(1,sigma,'.','MarkerSize',20, "MarkerEdgeColor","#4DBEEE", ...
    "MarkerFaceColor","#4DBEEE")
xlabel("$\mu$",'interpreter','latex')
ylabel("$\sigma^2$",'interpreter','latex')
ylim([3.8,4.2])
xlim([0.9,2])
Lgnd = legend("","$\mu$ = 0.1","", "$\mu$ = 1.0","","$\mu$ = 4.0","","$\mu$ = 10.0","True",'interpreter','latex');