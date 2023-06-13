
%% MU = 3, fixed step size
h = 0.1;
t0 = 0;
mu = 3;
tend = 4*mu;
x0 = [1.0; 1.0];

% using the Implicit euler method
options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @VanPolJac);
[X1,T1,function_calls1,hs1,rs1] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler",options);
function_calls111 = function_calls1;
hs111 = hs1;
rs111 = rs1;

h = 0.01;
% using the Implicit euler method
options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @VanPolJac);
[X2,T2,function_calls2,hs2,rs2] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler",options);
function_calls22 = function_calls2;
hs22 = hs2;
rs22 = rs2;

h = 0.001;
% using the Implicit euler method
options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @VanPolJac);
[X3,T3,function_calls3,hs3,rs3] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler",options);
function_calls33 = function_calls3;
hs33 = hs3;
rs33 = rs3;

options = odeset('RelTol',0.000001,'AbsTol',0.000001);
[Tode45,Xode45] = ode15s(@VanPol,[t0 tend],x0,options,mu);

[Tode15s,Xode15s]=ode45(@VanPol,[t0 tend],x0,options,mu);

%
fig = figure;

subplot(3, 3, 1)
plot(T1,X1(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('h = 0.1','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 2)
plot(T2,X2(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('h = 0.01','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 3)
plot(T3,X3(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('h = 0.001','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 4)
plot(T1,X1(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 5)
plot(T2,X2(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 6)
plot(T3,X3(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 7)
plot(X1(:,1),X1(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 8)
plot(X2(:,1),X2(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 9)
plot(X3(:,1),X3(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')


sgtitle('Implicit Euler with fixed step size, $\mu = 3$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = -0.02;
legend("Implicit Euler","ode45","ode15s",'NumColumns',3,'FontSize',16,'Interpreter','latex');
%% MU = 20, fixed step size
h = 0.1;
t0 = 0;
mu = 20;
tend = 4*mu;
x0 = [1.0; 1.0];

% using the Implicit euler method
options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @VanPolJac);
[X1,T1,function_calls4,hs1,rs4] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler",options);
hs4 = hs1;
h = 0.01;
% using the Implicit euler method
options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @VanPolJac);
[X2,T2,function_calls5,hs2,rs5] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler",options);
hs5 = hs2;

h = 0.001;
% using the Implicit euler method
options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @VanPolJac);
[X3,T3,function_calls6,hs3,rs6] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler",options);
hs6 = hs3;

options = odeset('RelTol',0.000001,'AbsTol',0.000001);
[Tode45,Xode45] = ode15s(@VanPol,[t0 tend],x0,options,mu);

[Tode15s,Xode15s]=ode45(@VanPol,[t0 tend],x0,options,mu);

%
fig = figure;

subplot(3, 3, 1)
plot(T1,X1(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('h = 0.1','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 2)
plot(T2,X2(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('h = 0.01','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 3)
plot(T3,X3(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('h = 0.001','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 4)
plot(T1,X1(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 5)
plot(T2,X2(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 6)
plot(T3,X3(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 7)
plot(X1(:,1),X1(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 8)
plot(X2(:,1),X2(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 9)
plot(X3(:,1),X3(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')


sgtitle('Implicit Euler with fixed step size, $\mu = 20$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = -0.02;
legend("Implicit Euler","ode45","ode15s",'NumColumns',3,'FontSize',16,'Interpreter','latex');

%% MU = 3, adaptive step size
h = 0.1;
t0 = 0;
mu = 3;
tend = 4*mu;
x0 = [1.0; 1.0];
Atol = 10^(-2);
Rtol = Atol;

% using the Implicit euler method
options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol, 'Jac' , @VanPolJac);
[X1,T1,function_calls7,hs1,rs7] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler", options);
hs7 = hs1;

Atol = 10^(-4);
Rtol = Atol;
% using the Implicit euler method
options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol, 'Jac' , @VanPolJac);
[X2,T2,function_calls8,hs2,rs8] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler", options);
hs8 = hs2;

Atol = 10^(-6);
Rtol = Atol;
% using the Implicit euler method
options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol, 'Jac' , @VanPolJac);
[X3,T3,function_calls9,hs3,rs9] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler", options);
hs9 = hs3;

options = odeset('RelTol',0.000001,'AbsTol',0.000001);
[Tode45,Xode45] = ode15s(@VanPol,[t0 tend],x0,options,mu);

[Tode15s,Xode15s]=ode45(@VanPol,[t0 tend],x0,options,mu);

%
fig = figure;

subplot(3, 3, 1)
plot(T1,X1(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-2}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 2)
plot(T2,X2(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-4}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 3)
plot(T3,X3(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-6}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 4)
plot(T1,X1(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 5)
plot(T2,X2(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 6)
plot(T3,X3(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 7)
plot(X1(:,1),X1(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 8)
plot(X2(:,1),X2(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 9)
plot(X3(:,1),X3(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')


sgtitle('Implicit Euler with adaptive step size, $\mu = 3$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = -0.02;
legend("Implicit Euler","ode45","ode15s",'NumColumns',3,'FontSize',16,'Interpreter','latex');


% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = -0.02;
legend("Implicit Euler","ode45","ode15s",'NumColumns',3,'FontSize',16,'Interpreter','latex');


fig = figure;
T1 = T1(2:end-1);
T2 = T2(2:end-1);
T3 = T3(2:end-1);
hs1 = hs1(:,1:end-1);
hs2 = hs2(:,1:end-1);
hs3 = hs3(:,1:end-1);
subplot(1, 3, 1)
semilogy(T1,hs1(1,:),'LineWidth', 1.5)
hold on
semilogy(T1(hs1(2,:)>1),hs1(1,hs1(2,:)>1),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-2}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

subplot(1, 3, 2)
semilogy(T2,hs2(1,:),'LineWidth', 1.5)
hold on
semilogy(T2(hs2(2,:)>1),hs2(1,hs2(2,:)>1),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-4}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

subplot(1, 3, 3)
semilogy(T3,hs3(1,:),'LineWidth', 1.5)
hold on
semilogy(T3(hs3(2,:)>1),hs3(1,hs3(2,:)>1),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-6}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

sgtitle('Implicit Euler with adaptive step size, $\mu = 3$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = 0;
legend("Step size","Rejections",'NumColumns',2,'FontSize',16,'Interpreter','latex');

%% MU = 20, adaptive step size
h = 0.1;
t0 = 0;
mu = 20;
tend = 4*mu;
x0 = [1.0; 1.0];
Atol = 10^(-2);
Rtol = Atol;

% using the Implicit euler method
options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol, 'Jac' , @VanPolJac);
[X1,T1,function_calls10,hs1,rs10] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler", options);
hs10 = hs1;

Atol = 10^(-4);
Rtol = Atol;
% using the Implicit euler method
options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol, 'Jac' , @VanPolJac);
[X2,T2,function_calls11,hs2,rs11] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler", options);
hs11 = hs2;

Atol = 10^(-6);
Rtol = Atol;
% using the Implicit euler method
options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol, 'Jac' , @VanPolJac);
[X3,T3,function_calls12,hs3,rs12] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler", options);
hs12 = hs3;

options = odeset('RelTol',0.000001,'AbsTol',0.000001);
[Tode45,Xode45] = ode15s(@VanPol,[t0 tend],x0,options,mu);

[Tode15s,Xode15s]=ode45(@VanPol,[t0 tend],x0,options,mu);

%
fig = figure;

subplot(3, 3, 1)
plot(T1,X1(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-2}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 2)
plot(T2,X2(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-4}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 3)
plot(T3,X3(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-6}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 4)
plot(T1,X1(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 5)
plot(T2,X2(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 6)
plot(T3,X3(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 7)
plot(X1(:,1),X1(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 8)
plot(X2(:,1),X2(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 9)
plot(X3(:,1),X3(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')


sgtitle('Implicit Euler with adaptive step size, $\mu = 20$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = -0.02;
legend("Implicit Euler","ode45","ode15s",'NumColumns',3,'FontSize',16,'Interpreter','latex');


fig = figure;
T1 = T1(2:end-1);
T2 = T2(2:end-1);
T3 = T3(2:end-1);
hs1 = hs1(:,1:end-1);
hs2 = hs2(:,1:end-1);
hs3 = hs3(:,1:end-1);
subplot(1, 3, 1)
semilogy(T1,hs1(1,:),'LineWidth', 1.5)
hold on
semilogy(T1(hs1(2,:)>1),hs1(1,hs1(2,:)>1),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-2}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

subplot(1, 3, 2)
semilogy(T2,hs2(1,:),'LineWidth', 1.5)
hold on
semilogy(T2(hs2(2,:)>1),hs2(1,hs2(2,:)>1),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-4}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

subplot(1, 3, 3)
semilogy(T3,hs3(1,:),'LineWidth', 1.5)
hold on
semilogy(T3(hs3(2,:)>1),hs3(1,hs3(2,:)>1),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-6}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

sgtitle('Implicit Euler with adaptive step size, $\mu = 20$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = 0;
legend("Step size","Rejections",'NumColumns',2,'FontSize',16,'Interpreter','latex');


%% Timing
h = 0.0001;
t0 = 0;
mu = 20;
tend = 4*mu;
x0 = [1.0; 1.0];

% using the Implicit euler method
options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @VanPolJac);
[X1,T1,function_calls66,hs66,rs66] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "Implicit Euler",options);

Stiff = zeros(3,7);
Stiff(1,1) = rs4;
Stiff(2,1) = function_calls4;
Stiff(3,1) = length(hs4);

Stiff(1,2) = rs5;
Stiff(2,2) = function_calls5;
Stiff(3,2) = length(hs5);

Stiff(1,3) = rs6;
Stiff(2,3) = function_calls6;
Stiff(3,3) = length(hs6);

Stiff(1,4) = rs66;
Stiff(2,4) = function_calls66;
Stiff(3,4) = length(hs66);

Stiff(1,5) = rs10;
Stiff(2,5) = function_calls10;
Stiff(3,5) = sum(hs10(2,:));

Stiff(1,6) = rs11;
Stiff(2,6) = function_calls11;
Stiff(3,6) = sum(hs11(2,:));

Stiff(1,7) = rs12;
Stiff(2,7) = function_calls12;
Stiff(3,7) = sum(hs12(2,:));

NonStiff = zeros(3,6);

NonStiff(1,1) = rs111;
NonStiff(2,1) = function_calls111;
NonStiff(3,1) = length(hs111);

NonStiff(1,2) = rs22;
NonStiff(2,2) = function_calls22;
NonStiff(3,2) = length(hs22);

NonStiff(1,3) = rs33;
NonStiff(2,3) = function_calls33;
NonStiff(3,3) = length(hs33);

NonStiff(1,4) = rs7;
NonStiff(2,4) = function_calls7;
NonStiff(3,4) = sum(hs7(2,:));

NonStiff(1,5) = rs8;
NonStiff(2,5) = function_calls8;
NonStiff(3,5) = sum(hs8(2,:));

NonStiff(1,6) = rs9;
NonStiff(2,6) = function_calls9;
NonStiff(3,6) = sum(hs9(2,:));



