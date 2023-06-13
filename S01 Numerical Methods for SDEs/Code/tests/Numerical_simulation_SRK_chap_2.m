


addpath(genpath("..\sde_solvers"))
% Duffing van der Pol

% Time-span
T = 2;
num = 3;
test = 14;
dt = 2^-num;
t = 0:dt:T;
dt_test = 2^(-test);
t_test = 0:dt_test:T;

% Brownian motion
rng(2);          % Set random seed
R_test_1 = randn(1,length(t_test));
R_test = cumsum(R_test_1)*sqrt(dt_test);
R = R_test(1:2^(test-num):end);
%R = [sum(reshape(R_test(1:end-round(2^(test-num)*(length(t_test)/2^(test-num)-floor(length(t_test)/2^(test-num))))),2^(test-num),floor(length(t_test)/2^(test-num)))),0];

figure; hold on 
set(gca,'TickLabelInterpreter','latex')
plot(t,R)
plot(t_test,R_test)
grid()
title("Brownian motion path",'interpreter','latex', FontSize=16)
legend('Low resolution','High resolution','interpreter','latex','Location','southeast')

% Parameters
alpha = 0.2;
beta = 1.5;

% The model
f = @(x,t) alpha*x;
g = @(x,t) beta*x;

x0 = 1/2;

% ODE
figure; hold on
set(gca,'TickLabelInterpreter','latex')

% Analytical solution
x = zeros(1,length(t_test));
x(1) = x0;
for i = 2:length(t_test)
    x(i) = x(i-1)*exp((alpha-beta^2/2)*dt_test+beta*(R_test(i)-R_test(i-1)));
end
plot(t_test,x);

plot(t,x(1:2^(test-num):end),'.','Color','green', 'MarkerSize',25);%[0.4660 0.6740 0.1880]



xlabel('$t$','interpreter','latex'); ylabel('$X_t$','interpreter','latex')
grid()

x_EM = eulermaruyama(f,g,t,x0,R);
plot(t(1:end),x_EM(1:end),'-o', 'MarkerSize',5,'LineWidth',1.1);

x_srkS10 = srkS10scalarnoise(f,g,t,x0,R);
plot(t(1:end),x_srkS10(1:end),'-o', 'MarkerSize',5,'LineWidth',1.1);

x_srkS15 = srkS15scalarnoise(f,g,t,x0,R);
plot(t(1:end),x_srkS15(1:end),'-o', 'Color',	"#D95319", 'MarkerSize',5,'LineWidth',1.1);

legend('Analytical path','Analytical points','EM','srkS10','srkS15','Location','northeast','interpreter','latex')
hold off



%%
rng(3);   
halfs = 7;
N = 1000;
errors_EM = zeros(N,halfs);
errors_srkS10 = zeros(N,halfs);
errors_srkS15 = zeros(N,halfs);


for k = 1:N
% Time-span
sum_x = 0;
sum_x_EM =0;
sum_x_srkS10 =0;
sum_x_srkS15 =0;

T = 5;
test = 14;
dt_test = 2^(-test);
t_test = 0:dt_test:T;
R_test_1 = randn(1,length(t_test));
R_test = cumsum(R_test_1)*sqrt(dt_test);

x = zeros(1,length(t_test));
x(1) = x0;
for i = 2:length(t_test)
    x(i) = x(i-1)*exp((alpha-beta^2/2)*dt_test+beta*(R_test(i)-R_test(i-1)));
end

for j = 1:halfs
    dt = 2^-(j);
    t = 0:dt:T;
    R = R_test(1:2^(test-j):end);
    % Brownian motion
    
    
    x_EM = eulermaruyama(f,g,t,x0,R);
    
    x_srkS10 = srkS10scalarnoise(f,g,t,x0,R);
        
    x_srkS15 = srkS15scalarnoise(f,g,t,x0,R);
    
    errors_srkS15(k,j) = abs(x_srkS15(end)-x(end));
    errors_srkS10(k,j) = abs(x_srkS10(end)-x(end));
    errors_EM(k,j) = abs(x_EM(end)-x(end));


end
    if mod(k,1000) == 0
        fprintf('At iteration %d...\n',k);
    end
end


%%
dts = 2.^(-halfs:-1);
p_srk15 = polyfit(log(flip(dts)),log(mean(errors_srkS15,1)),1);
p_srk1 = polyfit(log(flip(dts)),log(mean(errors_srkS10,1)),1);
p_em = polyfit(log(flip(dts)),log(mean(errors_EM,1)),1);

figure; hold on
set(gca,'TickLabelInterpreter','latex')

scatter(log(flip(dts)),log(mean(errors_EM,1)),"Color",	"#0072BD")
plot(log(flip(dts)),polyval(p_em,log(flip(dts))),"Color",	"#0072BD")

scatter(log(flip(dts)),log(mean(errors_srkS10,1)),"Color",	"#EDB120")
plot(log(flip(dts)),polyval(p_srk1,log(flip(dts))),"Color",	"#EDB120")

scatter(log(flip(dts)),log(mean(errors_srkS15,1)),"Color",	"#77AC30")
plot(log(flip(dts)),polyval(p_srk15,log(flip(dts))),"Color",	"#77AC30")

legend('EM','','srkS10','','srkS15','','Location','southeast','interpreter','latex')
xlabel('log$(dt)$','interpreter','latex'); ylabel('log(Error)','interpreter','latex')
title("Strong order", 'FontSize',16,'interpreter','latex')
grid on




p_srk15 = polyfit(log(flip(dts)),log(mean(errors_srkS15,1)),1);
p_srk1 = polyfit(log(flip(dts)),log(mean(errors_srkS10,1)),1);
p_em = polyfit(log(flip(dts)),log(mean(errors_EM,1)),1);

X = ['Observed strong order of Euler Maruyama is ',num2str(p_em(1))];
disp(X)

X = ['Observed strong order of SRK 1.0 is ',num2str(p_srk1(1))];
disp(X)

X = ['Observed strong order of SRK 1.5 is ',num2str(p_srk15(1))];
disp(X)

