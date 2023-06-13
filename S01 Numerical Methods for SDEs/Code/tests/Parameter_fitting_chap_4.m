addpath(genpath("..\sde_solvers"))
addpath(genpath("..\data"))
close all

T = readtable('data_for_parm_estimation_short_series.csv');
Futures = readtable('futures_data.csv');
Data = T;

M1vec = zeros(10000,1);
M1data = string(zeros(10000,1));

idx = zeros(10000,1);

j = 1;
for i=1:length(Futures{:,7})
    if strcmp(Futures{i,7}{1},'M1')
        M1data(i) = Futures{i,1};
        M1vec(i) = Futures{i,6};
        idx(j) = i;
        j = j+1;
    end

end

len = length(Data{:,1});
len = len -283;
idx = idx(1:len);



figure
plot(Futures{idx,1},Futures{idx,6})
grid()
ylabel("EUR/MWh",'interpreter','latex')
xlabel("Date",'interpreter','latex')

figure
plot(Data{:,1},Data{:,2})
grid()
ylabel("EUR/MWh",'interpreter','latex')
xlabel("Date",'interpreter','latex')

%

%OU_data = T{:,4}-T{:,2}-T{:,3}-T{:,5}-T{:,6}-T{:,7};
OU_data = T{:,7};

m_data = mean(OU_data);
var_data = var(OU_data);

norm_back = @(x) (x+m_data)*var_data;

%OU_data = (OU_data-m_data)/var_data;

dt = 1;

% Parameter estimation
% SRA3 MLE
a_f_add = @(r,sigma) (-dt^3*r^3+3*dt^2*r^2-6*dt*r+6)/6;
s_f_add = @(r,sigma) sigma*dt*(dt^2*r^2-3*dt*r+3)/3;

ell_SRA3 = @(p) 0.5 * sum(log(2*pi*s_f_add(p(1),p(2))) + 1/s_f_add(p(1),p(2)) * (OU_data(2:end) - a_f_add(p(1),p(2)) * OU_data(1:end-1)).^2);

p_SRA3 = fminsearch(ell_SRA3,[2 2]);

[p_SRA3_2,fval,exitflag,output,grad,hessian] = fminunc(ell_SRA3,[p_SRA3(1) p_SRA3(2)]);

var_SRA3 = inv(hessian);


% Exact MLE
a_fe = @(r,sigma) exp(-r * dt);
s_fe = @(r,sigma) sigma/(2*r) * (1 - exp(-2*r*dt));

ell = @(p) 0.5 * sum(log(2*pi*s_fe(p(1),p(2))) + 1/s_fe(p(1),p(2)) * (OU_data(2:end) - a_fe(p(1),p(2)) * OU_data(1:end-1)).^2);

p_e = fminsearch(ell,[2 2]);


% Euler-Maruyama MLE
a_f = @(r,sigma) (1 - 1*r * dt);
s_f = @(r,sigma) sigma*dt;

ell_em = @(p) 0.5 * sum(log(2*pi*s_f(p(1),p(2))) + 1/s_f(p(1),p(2)) * (OU_data(2:end) - a_f(p(1),p(2)) * OU_data(1:end-1)).^2);

p_em = fminsearch(ell_em,[2 2]);

disp('#### exact #####')
disp(p_e(1)/(1/dt))
disp(sqrt(p_e(2)/(1/dt)))

disp('#### EM #####')
disp(p_em(1)/(1/dt))
disp(sqrt(p_em(2)/(1/dt)))

disp('#### SRA3 #####')
disp(p_SRA3_2(1)/(1/dt))
disp(sqrt(p_SRA3_2(2)/(1/dt)))

disp("SRA3 from fminunc")
disp(p_SRA3(1)/(1/dt))
disp(sqrt(p_SRA3(2)/(1/dt)))

disp("Conf. int. for SRA3")
disp([p_SRA3(1)-1.96*sqrt(var_SRA3(1,1)),p_SRA3(1),p_SRA3(1)+1.96*sqrt(var_SRA3(1,1))])
disp([sqrt(p_SRA3(2)-1.96*sqrt(var_SRA3(2,2))),sqrt(p_SRA3(2)),sqrt(p_SRA3(2)+1.96*sqrt(var_SRA3(2,2)))])

% Reconstruction

r_sra = p_SRA3(1)/(1/dt);
sigma_sra = p_SRA3(2)/(1/dt);

r = p_e(1)/(1/dt);
sigma = p_e(2)/(1/dt);

steps = length(OU_data);
dt = 1;

X_sra = zeros(1,steps);
T_sra = zeros(1,steps);

X = zeros(1,steps);
T = zeros(1,steps);

x = OU_data(1);
x_sra = OU_data(1);
t = 0;

A_sra = exp(-r_sra*dt);
Sigma_sra = (sigma_sra - sigma_sra*exp(-2*dt*r_sra))/(2*r_sra);

r_sim = r;
A = exp(-r_sim*dt);
Sigma = (sigma - sigma*exp(-2*dt*r_sim))/(2*r_sim);
for k=1:steps
    X_sra(k) = x_sra;
    T_sra(k) = t;
    rnd = randn;
    x_sra = A_sra*x_sra + sqrt(Sigma_sra) * rnd;

    X(k) = x;
    T(k) = t;
    x = A*x + sqrt(Sigma) * rnd;
    t = t + dt;
end

%
figure
plot(Data{:,1},OU_data)
grid()
ylim([-20,30])
ylabel("EUR/MWh",'interpreter','latex')
title("Original",'interpreter','latex')


figure
plot(Data{:,1},X_sra)
grid()
ylim([-20,30])
ylabel("EUR/MWh",'interpreter','latex')
title("SRA3",'interpreter','latex')

figure
plot(Data{:,1},X_sra)
hold on
plot(Data{:,1},OU_data)
grid()
ylabel("EUR/MWh",'interpreter','latex')
title("Compare",'interpreter','latex')
legend("Recreation","Original",'interpreter','latex');

