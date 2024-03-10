clear
rng(12,'twister')
%%
%Number of steps
Nsteps = 100;
Nmc=1000; %Number of Monte Carlo runs
%f=1.57542e9/(2*pi); % 1575.42MHz (f的单位是rad)
f=2;
N_it=2;
N_x=4;
BW=20*pi;% 10Hz
load('scintDat.mat');
magnitude_full = abs(zkhist);
magnitude = magnitude_full (1:Nsteps);
phase_full = angle(zkhist);
phase = phase_full (1:Nsteps);
x_ini=[pi/2, f, 0.5, phase(1)]';
P=diag([(pi^2)/3 1 1 0.08]); % 0.08 is from Valis paper
T=0.01;
% T=1/(5*f);%sampling period
% G = [T^3/3; T^2/2; T];
% G_matrix = G*G.' * 1;
% Q = [G_matrix, zeros(3, 1); zeros(1, 3), 0.08];%covariance matrix
Q=0.1*diag([0.1 0.1 0.1 0.1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);

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
[y,y_3,RMSE_2,RMSE_3] = PLL_RMSE(BW,f,Nsteps,T,y_measure_re,x_truth);
%% Draw
draw_filtered(Nsteps,x_u_series_EKF,x_u_series_UKF,x_u_series_IPLF,x_truth,y_measure_re,y,y_3)
%% MC simulation
independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol_EKF=zeros(Nsteps,Nmc);
RMSE_tol_UKF=zeros(Nsteps,Nmc);
RMSE_tol_IPLF=zeros(Nsteps,Nmc);
RMSE_tot_2=zeros(Nsteps,Nmc);
RMSE_tot_3=zeros(Nsteps,Nmc);
% RMSE_tol_IPLF=zeros(Nsteps,Nmc);
for i=1:Nmc      
    [y_measure_re,y_theta,x_truth_phase,y_measure_mc,x_truth] = generate_truth_PLL(Nsteps,x_ini,Q,R,T);
    %EKF
    [x_u_series,RMSE_tol_EKF(:,i)] = EKF(Nsteps,x_ini,P,R,Q,F,y_measure_mc,x_truth, phase,magnitude);
    %UKF
    [x_u_series_UKF,RMSE_tol_UKF(:,i)] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure_mc,N_x,x_truth,magnitude);
    %IPLF
    [x_u_series_IPLF,RMSE_tol_IPLF(:,i)] = IPLF(Nsteps,x_ini,P,R,Q,F,N_x,x_truth,N_it,y_measure_mc);
    %PLL
    [y,y_3,RMSE_tot_2(:,i),RMSE_tot_3(:,i)] = PLL_RMSE(BW,f,Nsteps,T,y_measure_re,x_truth);

end

rmse_EKF=sum(RMSE_tol_EKF,2)/Nmc;
rmse_UKF=sum(RMSE_tol_UKF,2)/Nmc;
rmse_IPLF=sum(RMSE_tol_IPLF,2)/Nmc;
rmse_2=sum(RMSE_tot_2,2)/Nmc;
rmse_3=sum(RMSE_tot_3,2)/Nmc;

%% Draw RMSE   
    figure(2)
    orange = [1, 0.5, 0];% UKF
    black = [0 0 0];% IPLF
    blue = [0.21 0.35 1];% EKF 
    yellow = [1, 0.84, 0];
    purple = [0.54, 0.17, 0.89];% 2 PLL
    white = [1 1 1];
    green = [0.1 0.8 0.5];% 3 PLL
    
    lightgrey = [0.94 0.94 0.94]; % color definition

    plot((1:Nsteps),rmse_EKF,'--','Color',blue,'LineWidth',2)
    hold on     
    plot((1:Nsteps),rmse_UKF,'--','Color',orange,'LineWidth',2)
    hold on
    plot((1:Nsteps),rmse_IPLF,'--','Color',black,'LineWidth',2)
    hold on
    plot((1:Nsteps),rmse_3,'.-','Color',green,'LineWidth',2)
    hold on
    plot((1:Nsteps),rmse_2,'.-','Color',purple,'LineWidth',2)

    h1 = legend('RMSE(EKF)','RMSE(UKF)','RMSE(IPLF)','RMSE(3rd-order-PLL)','RMSE(2nd-order-PLL)','Location','northeast'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',white) %filling colour of legend
    set(h1,'Box','on') %Remove outer frame of legend
    grid on
    box off
    title('RMSE comparison')
    ylabel('RMS phase error [rad]')
    xlabel('Time step')
    grid on
axis([ 0 Nsteps 0 max(rmse_2)+0.05]) 

%% Output mean RMSE
averageValue_UKF = mean(rmse_UKF(:));
fprintf('RMSE_UKF=%0.5f\n', averageValue_UKF);
averageValue_EKF = mean(rmse_EKF(:));
fprintf('RMSE_EKF=%0.5f\n', averageValue_EKF);
averageValue_IPLF = mean(rmse_IPLF(:));
fprintf('RMSE_IPLF=%0.5f\n', averageValue_IPLF);
averageValue_PLL_2 = mean(rmse_2);
fprintf('RMSE_2_PLL=%0.5f\n', averageValue_PLL_2);
averageValue_PLL_3 = mean(rmse_3);
fprintf('RMSE_3_PLL=%0.5f\n', averageValue_PLL_3);




