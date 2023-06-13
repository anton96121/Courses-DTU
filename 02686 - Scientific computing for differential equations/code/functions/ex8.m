h = 0.1;
t0 = 0;
mu = 5:5:120;
x0 = [1.0; 1.0];
Atol = 10^(-6);
Rtol = Atol;
tests = length(mu);

function_calls_ESDIRK = zeros(tests,1);
jacobian_calls_ESDIRK = zeros(tests,1);
steps_ESDIRK = zeros(tests,1);
rs_ESDIRK = zeros(tests,1);
reject_ESDIRK = zeros(tests,1);

function_calls_DOPRI = zeros(tests,1);
rs_DOPRI = zeros(tests,1);
steps_DOPRI = zeros(tests,1);
reject_DOPRI = zeros(tests,1);

for j = 1:10
for i = 1:tests
    mu_i = mu(i);
    tend = 20*mu_i;
    % using the ESDIRK23method
    options = struct('Jac' , @VanPolJac,'Rtol',Rtol,'Atol',Atol);
    [X1,T1,function_calls1,hs1,rs1] = ODEsolver(@VanPol,[mu_i],h,t0,tend,x0, "ESDIRK23", options);
    function_calls_ESDIRK(i) =  function_calls1.nFun;
    jacobian_calls_ESDIRK(i) = function_calls1.nJac;
    rs_ESDIRK(i) = rs_ESDIRK(i)+rs1;
    steps_ESDIRK(i) = function_calls1.nStep;
    reject_ESDIRK(i) = function_calls1.nFail;
    % using the Dormand-Prince 5(4)method
    options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol);
    [X2,T2,function_calls2,h21,rs2] = ODEsolver(@VanPol,[mu_i],h,t0,tend,x0, "DOPRI54", options);
    function_calls_DOPRI(i) = function_calls2;
    rs_DOPRI(i) =rs_DOPRI(i)+ rs2;
    steps_DOPRI(i) = sum(h21(2,:));
    reject_DOPRI(i) = sum(h21(2,:)-1);
    disp(i)
end
disp("j runde")
disp(j)
disp("####")
end

%%

figure
subplot(2, 2, 1)
plot(mu,rs_DOPRI/10)
hold on 
plot(mu,rs_ESDIRK/10)
xlabel("$\mu$",'FontSize',16,'Interpreter','latex')
ylabel("Time [s]",'FontSize',16,'Interpreter','latex')
legend("Dormand-Prince 5(4)","ESDIRK23",'FontSize',12,'Interpreter','latex','Location','northwest')
title("Compute time",'FontSize',20,'Interpreter','latex')

subplot(2, 2, 2)
semilogy(mu,steps_DOPRI)
hold on 
semilogy(mu,steps_ESDIRK)
xlabel("$\mu$",'FontSize',16,'Interpreter','latex')
ylabel("Steps",'FontSize',16,'Interpreter','latex')
legend("Dormand-Prince 5(4)","ESDIRK23",'FontSize',12,'Interpreter','latex','Location','northwest')
title("Steps used",'FontSize',20,'Interpreter','latex')

subplot(2, 2, 3)
plot(mu,reject_DOPRI)
hold on 
plot(mu,reject_ESDIRK)
xlabel("$\mu$",'FontSize',16,'Interpreter','latex')
ylabel("Rejected steps",'FontSize',16,'Interpreter','latex')
legend("Dormand-Prince 5(4)","ESDIRK23",'FontSize',12,'Interpreter','latex','Location','northwest')
title("Rejected steps",'FontSize',20,'Interpreter','latex')



subplot(2, 2, 4)
semilogy(mu,function_calls_DOPRI)
hold on 
semilogy(mu,function_calls_ESDIRK)
semilogy(mu,jacobian_calls_ESDIRK)
semilogy(mu,function_calls_ESDIRK+jacobian_calls_ESDIRK)
xlabel("$\mu$",'FontSize',16,'Interpreter','latex')
ylabel("Function calls",'FontSize',16,'Interpreter','latex')
legend("Dormand-Prince 5(4)","ESDIRK23, gradient calls","ESDIRK23, jacobian calls", "ESDIRK23, total",'FontSize',12,'Interpreter','latex','Location','northwest')
title("Function calls",'FontSize',20,'Interpreter','latex')

%%
h = 0.1;
t0 = 0;
mu = 120;
x0 = [1.0; 1.0];
Atol = 10^(-6);
Rtol = Atol;

% using the ESDIRK23method
options = struct('Jac' , @VanPolJac,'Rtol',Rtol,'Atol',Atol);
[X1,T1,function_calls1,hs1,rs1] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "ESDIRK23", options);

% using the Dormand-Prince 5(4)method
options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol);
[X2,T2,function_calls2,h21,rs2] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "DOPRI54", options);

options = odeset('RelTol',0.000001,'AbsTol',0.000001);
[Tode45,Xode45] = ode15s(@VanPol,[t0 tend],x0,options,mu);

figure
plot(T1,X1(:,1),'LineWidth', 1.5)
hold on 
plot(T2,X2(:,1),'LineWidth', 1.5)
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
title('Atol = Rtol = $10^{-6}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')
legend("ESDIRK23", "DOPRI54","ode15s")

%%
h = 0.1;
t0 = 0;
mu = 2:2:20;
x0 = [1.0; 1.0];
Atol = 10^(-6);
Rtol = Atol;
tests = length(mu);

function_calls_RK4 = zeros(tests,1);
rs_RK4 = zeros(tests,1);
steps_RK4 = zeros(tests,1);
reject_RK4 = zeros(tests,1);

function_calls_DOPRI = zeros(tests,1);
rs_DOPRI = zeros(tests,1);
steps_DOPRI = zeros(tests,1);
reject_DOPRI = zeros(tests,1);

for j = 1:10
for i = 1:tests
    mu_i = mu(i);
    tend = 20*mu_i;
    % using the ESDIRK23method
    options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol);
    [X2,T2,function_calls2,h21,rs2] = ODEsolver(@VanPol,[mu_i],h,t0,tend,x0, "RK4", options);
    function_calls_RK4(i) = function_calls2;
    rs_RK4(i) =rs_DOPRI(i)+ rs2;
    steps_RK4(i) = sum(h21(2,:));
    reject_RK4(i) = sum(h21(2,:)-1);
    disp(i)
    % using the Dormand-Prince 5(4)method
    options = struct('step_control',true, 'initialStepSize',true,'control_type', "PID",'Rtol',Rtol,'Atol',Atol);
    [X2,T2,function_calls2,h21,rs2] = ODEsolver(@VanPol,[mu_i],h,t0,tend,x0, "DOPRI54", options);
    function_calls_DOPRI(i) = function_calls2;
    rs_DOPRI(i) =rs_DOPRI(i)+ rs2;
    steps_DOPRI(i) = sum(h21(2,:));
    reject_DOPRI(i) = sum(h21(2,:)-1);
    disp(i)
end
disp("j runde")
disp(j)
disp("####")
end

%%
figure
subplot(2,2,1)
plot(mu,rs_DOPRI/10)
hold on 
plot(mu,rs_RK4/10)
xlabel("$\mu$",'FontSize',16,'Interpreter','latex')
ylabel("Time [s]",'FontSize',16,'Interpreter','latex')
legend("Dormand-Prince 5(4)","RK4",'FontSize',12,'Interpreter','latex','Location','northwest')
title("Compute time",'FontSize',20,'Interpreter','latex')

subplot(2,2,2)
plot(mu,function_calls_DOPRI)
hold on 
plot(mu,function_calls_RK4)
xlabel("$\mu$",'FontSize',16,'Interpreter','latex')
ylabel("Function calls",'FontSize',16,'Interpreter','latex')
legend("Dormand-Prince 5(4)","RK4",'FontSize',12,'Interpreter','latex','Location','northwest')
title("Function calls",'FontSize',20,'Interpreter','latex')

subplot(2,2,3)
plot(mu,steps_DOPRI)
hold on 
plot(mu,steps_RK4)
xlabel("$\mu$",'FontSize',16,'Interpreter','latex')
ylabel("Steps",'FontSize',16,'Interpreter','latex')
legend("Dormand-Prince 5(4)","RK4",'FontSize',12,'Interpreter','latex','Location','northwest')
title("Steps",'FontSize',20,'Interpreter','latex')
