clear
rng(12,'twister')
%%
%Number of steps
Nsteps = 100;
Nmc=1000; %Number of Monte Carlo runs
f=1.57542e9/(2*pi); % 1575.42MHz (f的单位是rad)
N_it=2;
N_x=4;
BW=40*pi;
load('scintDat.mat');
magnitude_full = abs(zkhist);
magnitude = magnitude_full (1:100);
phase_full = angle(zkhist);
phase = phase_full (1:100);
x_ini=[pi/2, f, 0, phase(1)]';
P=diag([(pi^2)/3 1 1 0.08]); % 0.08 is from Valis paper
Q=0.1*diag([0.1 0.1 0.1 1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);


T=1/(5*f);%sampling period
beta = 0.914;
F=[1 T (T^2)/2 0; 0 1 T 0; 0 0 T 0; 0 0 0 beta];

[y_measure_re,y_theta,x_truth_phase,y_measure,x_truth] = generate_truth_PLL(Nsteps,x_ini,Q,R,T);

%% EKF-AR
[x_u_series_EKF,RMSE_EKF] = EKF(Nsteps,x_ini,P,R,Q,F,y_measure,x_truth,phase,magnitude);
%% UKF-AR
[x_u_series_UKF,RMSE_UKF] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure,N_x,x_truth,magnitude);
%% IPLF-AR
[x_u_series_IPLF,RMSE_IPLF] = IPLF(Nsteps,x_ini,P,R,Q,F,N_x,x_truth,N_it,y_measure);
%% PLL
% in <draw_filter.m>
%% Draw
[average_RMSE_2,average_RMSE_3,RMSE_2,RMSE_3]=draw_filtered(Nsteps,x_u_series_EKF,x_u_series_UKF,x_u_series_IPLF,x_truth,y_measure_re,BW,f);
%% MC simulation
independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol_EKF=zeros(Nsteps,Nmc);
RMSE_tol_UKF=zeros(Nsteps,Nmc);
RMSE_tol_IPLF=zeros(Nsteps,Nmc);
% RMSE_tol_IPLF=zeros(Nsteps,Nmc);
for i=1:Nmc      
    [y_measure_re,y_theta,x_truth_phase,y_measure_mc,x_truth] = generate_truth_PLL(Nsteps,x_ini,Q,R,T);
    %EKF
    [x_u_series,RMSE_tol_EKF(:,i)] = EKF(Nsteps,x_ini,P,R,Q,F,y_measure_mc,x_truth, phase,magnitude);
    %UKF
    [x_u_series_UKF,RMSE_tol_UKF(:,i)] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure_mc,N_x,x_truth,magnitude);
    %IPLF
    [x_u_series_IPLF,RMSE_tol_IPLF(:,i)] = IPLF(Nsteps,x_ini,P,R,Q,F,N_x,x_truth,N_it,y_measure_mc);

end

rmse_EKF=sum(RMSE_tol_EKF,2)/Nmc;
rmse_UKF=sum(RMSE_tol_UKF,2)/Nmc;
rmse_IPLF=sum(RMSE_tol_IPLF,2)/Nmc;
%% Draw RMSE   
    figure(2)
    orange = [1, 0.5, 0];% UKF
    black = [0 0 0];% IPLF
    blue = [0.21 0.35 1];% EKF 
    yellow = [1, 0.84, 0];% 3 PLL
    purple = [0.54, 0.17, 0.89];% 2 PLL
    white = [1 1 1];
    
    lightgrey = [0.94 0.94 0.94]; % color definition

    plot((1:Nsteps),rmse_EKF,'--','Color',blue,'LineWidth',2)
    hold on     
    plot((1:Nsteps),rmse_UKF,'--','Color',orange,'LineWidth',2)
    hold on
    plot((1:Nsteps),rmse_IPLF,'--','Color',black,'LineWidth',2)
    hold on
    plot((1:Nsteps),RMSE_3,'.-','Color',yellow,'LineWidth',2)
    hold on
    plot((1:Nsteps),RMSE_2,'.-','Color',purple,'LineWidth',2)

    h1 = legend('rmse(EKF)','rmse(UKF)','rmse(IPLF)','rmse(3rd-order-PLL)','rmse(2nd-order-PLL)','Location','southwest'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle',':');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',white) %filling colour of legend
    set(h1,'Box','on') %Remove outer frame of legend
    grid on
    box off
    title('RMSE comparison')
    ylabel('RMS phase error [rad]')
    xlabel('Time step')
    grid on
axis([ 0 Nsteps 0 0.12]) 
% axis([ 0 Nsteps 0 max(rmse_EKF)+0.05]) 
averageValue_EKF = mean(rmse_EKF(:))
averageValue_UKF = mean(rmse_UKF(:))
averageValue_IPLF = mean(rmse_IPLF(:))
averageValue_PLL_2 = mean(average_RMSE_2)
averageValue_PLL_3 = mean(average_RMSE_3)




