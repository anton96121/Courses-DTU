addpath(genpath("..\data"))
addpath(genpath("..\figures_report"))


parm_estim_data = readtable("../data/data_for_parm_estimation_short_series.csv");
%all_futures_data = readtable("..\data\futures_data.csv");

m1_future = readtable("../data/data_future_found/short_period/data_m1_short_period.csv");
m2_future = readtable("../data/data_future_found/short_period/data_m2_short_period.csv");
m3_future = readtable("../data/data_future_found/short_period/data_m3_short_period.csv");
m4_future = readtable("../data/data_future_found/short_period/data_m4_short_period.csv");
%idx_end = find(eex_data.date== '2009-01-01')-1;
%data = table2array(eex_data(1:idx_end,"price"));
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
%%

dt = 2^-5;
t_test = 0:dt:100;
rng(200);    

brownian = cumsum(randn(1,length(t_test)));

x_brow = zeros(1,length(t_test));
x_parm_change = zeros(1,length(t_test));
x_brow(1) = 0;
for i = 2:length(t_test)
    x_brow(i) = x_brow(i-1) + 1*(brownian(i)-brownian(i-1));
    x_parm_change(i) = x_parm_change(i-1) + 3*dt + 1*(brownian(i)-brownian(i-1));
    %x_brow(i) = x(i-1)*exp((alpha-beta^2/2)*dt_test+beta*(brownian(i)-brownian(i-1)));
end
figure()
hold on 
plot(t_test,x_brow);
plot(t_test,x_parm_change);
hold off
xlabel('time','interpreter','latex')
legend("$B_t$","$\bar{B}_t$", 'Interpreter','latex','location','northwest')
saveas(gcf, "../figures_report/application/figure_theory/brownian_and_scaled_brownian",'epsc');


%%
head(parm_estim_data,3)

%% Spot Price
figure()
hold on
plot(parm_estim_data{:,1},parm_estim_data{:,4})
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
hold off
saveas(gcf, "../figures_report/application/initial_overview/spot_price_entire_period",'epsc');
saveas(gcf, "../figures_report/application/initial_overview/spot_price_entire_period",'png');

%% Future for each M1, M2, M3, and M4
list_ms = ["M1" "M2" "M3" "M4"];
for i=1:4
    mask_mi =  all_futures_data.forward_month == list_ms(i);
    plot(all_futures_data{mask_mi,'notationtime'},all_futures_data{mask_mi,'price'})
    xlabel('date','interpreter','latex')
    ylabel('EUR/MWh','interpreter','latex')
    lable_image = strcat("../figures_report/application/initial_overview/future_price_entire_period_",list_ms(i));
    saveas(gcf, lable_image,'epsc');
    saveas(gcf, lable_image,'png');
end


%% Plot ACF, PACF and Variogram

%% Only the first lags

figure()
autocorr(parm_estim_data{:,4},NumLags=30)
xlabel('Lag','interpreter','latex')
ylabel('Sample Autocorrelation','interpreter','latex')
title('')
%title('Sample Autocorrelation Function','interpreter','latex')
saveas(gcf, "../figures_report/application/initial_overview/ACF_entire_period_30_lags",'epsc');
saveas(gcf, "../figures_report/application/initial_overview/ACF_entire_period_30_lags",'png');


%%
figure()
parcorr(parm_estim_data{:,4},NumLags=30)
xlabel('Lag','interpreter','latex')
ylabel('Sample Partial Autocorrelation','interpreter','latex')
title('')
%title('Sample Autocorrelation Function','interpreter','latex')
saveas(gcf, "../figures_report/application/initial_overview/PACF_entire_period_30_lags",'epsc');
saveas(gcf, "../figures_report/application/initial_overview/PACF_entire_period_30_lags",'png');


%%


figure()
autocorr(parm_estim_data{:,4},NumLags=600)
xlabel('Lag','interpreter','latex')
ylabel('Sample Autocorrelation','interpreter','latex')
title('')
%title('Sample Autocorrelation Function','interpreter','latex')
saveas(gcf, "../figures_report/application/initial_overview/ACF_entire_period",'epsc');
saveas(gcf, "../figures_report/application/initial_overview/ACF_entire_period",'png');

%% 
figure()
parcorr(parm_estim_data{:,4},NumLags=600)
xlabel('Lag','interpreter','latex')
ylabel('Sample Partial Autocorrelation','interpreter','latex')
title('')
%title('Sample Autocorrelation Function','interpreter','latex')
saveas(gcf, "../figures_report/application/initial_overview/PACF_entire_period",'epsc');
saveas(gcf, "../figures_report/application/initial_overview/PACF_entire_period",'png');

%%
[pxx,f] = periodogram(parm_estim_data{:,4},[],[],365);
figure()  
semilogy(f,pxx)
xlabel('Cycles/Year','interpreter','latex')
ylabel('Power','interpreter','latex')
%title('Periodogram of Spot Price','interpreter','latex')
saveas(gcf, "../figures_report/application/initial_overview/periodogram_entire_period",'epsc');
saveas(gcf, "../figures_report/application/initial_overview/periodogram_entire_period",'png');
%%
%{weekd=weekday(parm_estim_data{:,1}),parm_estim_data{:,4}}
[no,day_name] = weekday(parm_estim_data{:,1});
boxplot(parm_estim_data{:,4},{day_name})
ylim([0,200])
xlabel('Day of Week','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
saveas(gcf, "../figures_report/application/initial_overview/boxplot_day_of_week_entire_period",'epsc');
saveas(gcf, "../figures_report/application/initial_overview/boxplot_day_of_week_entire_period",'png');


%% deterministic function

determ_function = sum(parm_estim_data{:,["seasonal_7","seasonal_182","seasonal_365","trend"]},2);

figure()
hold on
plot(parm_estim_data{:,1},parm_estim_data{:,4})
plot(parm_estim_data{:,"date"},determ_function)
hold off 
ylim([0,150])
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
legend('S(t)','$\Lambda(t)$', 'Position',[0.35 0.75 0.1 0.1],'interpreter','latex')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/deterministic_function",'epsc')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/deterministic_function",'png')

%% deterministic function componentwise

figure()
plot(parm_estim_data{:,1},parm_estim_data{:,"trend"})
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/trend_component",'epsc')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/trend_component",'png')
%% Deterministic weekly profile

mask_week = 7:13;
week_days = categorical(["Mon", "Tue", "Wed", "Thur", "Fri", "Sat", "Sun"]);
figure()
plot(parm_estim_data{mask_week,"seasonal_7"})
xlabel('day of week','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
set(gca,'xticklabel',week_days.')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/weekly_component",'epsc')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/weekly_component",'png')

%% Deterministic semi-annual
idx_end = find(parm_estim_data{:,1}== '2002-08-01')-1;
mask_semi_annual = 1:idx_end;
figure()
plot(parm_estim_data{mask_semi_annual,"date"}, parm_estim_data{mask_semi_annual,"seasonal_182"})
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/semi-annual_component",'epsc')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/semi-annual_component",'png')

%% Deterministic annual
idx_end = find(parm_estim_data{:,1}== '2003-01-01')-1;
mask_annual = 1:idx_end;
figure()
plot(parm_estim_data{mask_annual,"date"}, parm_estim_data{mask_annual,"seasonal_365"})
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/annual_component",'epsc')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/annual_component",'png')

%% spike component

figure()
hold on
plot(parm_estim_data{:,1},parm_estim_data{:,'spike_vec'})
%plot(parm_estim_data{:,"date"},determ_function)
hold off 
%legend('S(t)','\Lambda(t)', 'Position',[0.35 0.75 0.1 0.1])
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/spike_signal_vector",'epsc')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/spike_signal_vector",'png')



figure()
hold on
plot(parm_estim_data{:,1},parm_estim_data{:,'init_data'}- parm_estim_data{:,'spike_vec'})
%plot(parm_estim_data{:,"date"},determ_function)
ylim([0,350])
hold off 
%legend('S(t)','\Lambda(t)', 'Position',[0.35 0.75 0.1 0.1])
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/spike_sub_init_data",'epsc')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/spike_sub_init_data",'png')


%% The Ornstein Uhlenbeck process in the Signal:

figure()
hold on
plot(parm_estim_data{:,1}, parm_estim_data{:,"data_estimate_brownian"})
%plot(parm_estim_data{:,"date"},determ_function)
hold off 
%legend('S(t)','\Lambda(t)', 'Position',[0.35 0.75 0.1 0.1])
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/ou_process_in_signal",'epsc')
saveas(gcf, "../figures_report/application/parts_of_fitted_model/ou_process_in_signal",'png')


%% M1 Analysis 

figure()
hold on
plot(m1_future{:,'notationtime'}, m1_future{:,"price"})
plot(m1_future{:,'notationtime'}, m1_future{:,"calc_under_q_theta"})
plot(m1_future{:,'notationtime'}, m1_future{:,"calc_under_p"})
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
legend("M1 - realized","M1 - under $\mathbf{Q_\theta}$","M1 - under $\mathbf{P}$", 'Interpreter','latex' ...
    ,'location','northwest')
grid()
hold off
saveas(gcf, "../figures_report/application/short_period_futures/m1_future_under_measures",'epsc')
saveas(gcf, "../figures_report/application/short_period_futures/m1_future_under_measures",'png')

%% M2 Analysis

figure()
hold on
plot(m2_future{:,'notationtime'}, m2_future{:,"price"})
plot(m2_future{:,'notationtime'}, m2_future{:,"calc_under_q_theta"})
plot(m2_future{:,'notationtime'}, m2_future{:,"calc_under_p"})
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
legend("M2 - realized","M2 - under $\mathbf{Q_\theta}$","M2 - under $\mathbf{P}$", 'Interpreter','latex' ...
    ,'location','northwest')
grid()
hold off

saveas(gcf, "../figures_report/application/short_period_futures/m2_future_under_measures",'epsc')
saveas(gcf, "../figures_report/application/short_period_futures/m2_future_under_measures",'png')

%% M3 

figure()
hold on
plot(m3_future{:,'notationtime'}, m3_future{:,"price"})
plot(m3_future{:,'notationtime'}, m3_future{:,"calc_under_q_theta"})
plot(m3_future{:,'notationtime'}, m3_future{:,"calc_under_p"})
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
legend("M3 - realized","M3 - under $\mathbf{Q_\theta}$","M3 - under $\mathbf{P}$", 'Interpreter','latex' ...
    ,'location','northwest')
grid()
hold off

saveas(gcf, "../figures_report/application/short_period_futures/m3_future_under_measures",'epsc')
saveas(gcf, "../figures_report/application/short_period_futures/m3_future_under_measures",'png')


%% M4

figure()
hold on
plot(m4_future{:,'notationtime'}, m4_future{:,"price"})
plot(m4_future{:,'notationtime'}, m4_future{:,"calc_under_q_theta"})
plot(m4_future{:,'notationtime'}, m4_future{:,"calc_under_p"})
xlabel('date','interpreter','latex')
ylabel('EUR/MWh','interpreter','latex')
legend("M4 - realized","M4 - under $\mathbf{Q_\theta}$","M4 - under $\mathbf{P}$", 'Interpreter','latex' ...
    ,'location','northwest')
grid()
hold off

saveas(gcf, "../figures_report/application/short_period_futures/m4_future_under_measures",'epsc')
saveas(gcf, "../figures_report/application/short_period_futures/m4_future_under_measures",'png')

%%