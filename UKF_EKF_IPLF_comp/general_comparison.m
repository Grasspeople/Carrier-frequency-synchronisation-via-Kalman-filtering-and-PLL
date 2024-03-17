clear
rng(12,'twister')
%% Parameters
Nsteps = 100; N_it = 5;                lambda=10^(-6)*3-3;
f=2;x_ini=[pi/2,f,0.5]';
P=diag([(pi^2)/3 1 1]);
Q=0.01*diag([0.01 0.1 0.1]);%covariance matrix %偷偷改的
R=diag([(pi/3)^2 (pi/3)^2]);
N_x=3;
BW=20*pi;%20Hz--40*pi///10Hz--20*pi///5Hz---10*pi

T=0.01;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];
Nmc=1000; %Number of Monte Carlo runs
independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
UKF_RMSE_tol=zeros(Nsteps,Nmc);
EKF_RMSE_tol=zeros(Nsteps,Nmc);
IPLF_RMSE_tol=zeros(Nsteps,Nmc);
RMSE_tot_2=zeros(Nsteps,Nmc);
RMSE_tot_3=zeros(Nsteps,Nmc);
for i=1:Nmc 
% Measurements
[y_measure_re,y_theta,x_truth_phase,y_measure,x_truth] = truth_comp(Nsteps,x_ini,Q,R,F);
[x_u_series_UKF,UKF_RMSE_tol(:,i)] = UKF_comp(Nsteps,x_ini,P,R,Q,F,y_measure,N_x,x_truth);
[x_u_series_EKF,EKF_RMSE_tol(:,i)] = EKF_comp(Nsteps,x_ini,P,R,Q,F,y_measure,x_truth);
[x_u_series_IPLF,IPLF_RMSE_tol(:,i)] = IPLF_comp(Nsteps,x_ini,P,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
[y,y_3,RMSE_tot_2(:,i),RMSE_tot_3(:,i)]=PLL_RMSE(Nsteps,f,BW,T,y_measure_re,x_truth_phase);
end
RMSE_UKF=sum(UKF_RMSE_tol,2)/Nmc;
RMSE_EKF=sum(EKF_RMSE_tol,2)/Nmc;
RMSE_IPLF=sum(IPLF_RMSE_tol,2)/Nmc;
RMSE_PLL_2=sum(RMSE_tot_2,2)/Nmc;
RMSE_PLL_3=sum(RMSE_tot_3,2)/Nmc;
%% Draw RMSE (whole)
    figure(1)
    orange = [1, 0.5, 0];% UKF
    black = [0 0 0];% IPLF
    blue = [0.21 0.35 1];% EKF 
    %yellow = [1, 0.84, 0];
    %yellow = [0.8, 0.7, 0.1];
    green = [0.1 0.8 0.5];% 3 PLL
    purple = [0.54, 0.17, 0.89];% 2 PLL
    white = [1 1 1];
    
    lightgrey = [0.94 0.94 0.94]; % color definition

    plot((1:Nsteps),RMSE_EKF,'--','Color',blue,'LineWidth',2)
    hold on     
    plot((1:Nsteps),RMSE_UKF,'--','Color',orange,'LineWidth',2)
    hold on
    plot((1:Nsteps),RMSE_IPLF,'--','Color',black,'LineWidth',2)
    hold on
    plot((1:Nsteps),RMSE_PLL_3,'.-','Color',green,'LineWidth',2)
    hold on
    plot((1:Nsteps),RMSE_PLL_2,'.-','Color',purple,'LineWidth',2)

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
    axis([ 0 Nsteps 0 max(RMSE_PLL_2)+0.05]) 
%      axis([ 0 Nsteps 0 0.3]) 
%% Output mean RMSE
averageValue_UKF = mean(RMSE_UKF(:));
fprintf('RMSE_UKF=%0.5f\n', averageValue_UKF);
averageValue_EKF = mean(RMSE_EKF(:));
fprintf('RMSE_EKF=%0.5f\n', averageValue_EKF);
averageValue_IPLF = mean(RMSE_IPLF(:));
fprintf('RMSE_IPLF=%0.5f\n', averageValue_IPLF);
averageValue_PLL_2 = mean(RMSE_PLL_2(:));
fprintf('RMSE_2_PLL=%0.5f\n', averageValue_PLL_2);
averageValue_PLL_3 = mean(RMSE_PLL_3(:));
fprintf('RMSE_3_PLL=%0.5f\n', averageValue_PLL_3);
%% Draw prediction
    orange = [1 0.34 0.20]; %UKF
    black = [0 0 0]; %IPLF
    blue = [0.21 0.35 1]; %EKF
    green = [0.1 0.8 0.5]; %truth
    purple = [0.75, 0.1, 0.75]; % 2 PLL
    yellow = [0.8, 0.7, 0.1];% 3 PLL
    brown = [0.65, 0.16, 0.16]; %observation
    lightgrey = [0.94 0.94 0.94]; %bg


    figure(2)
    plot((1:Nsteps),x_truth(1,:),'.-','Color',green,'LineWidth',3)
    hold on    
    plot((1:Nsteps),y_measure_re,'--','Color',brown,'LineWidth',2)
    hold on 
    plot((1:Nsteps),x_u_series_EKF(1,:),'.-','Color',blue,'LineWidth',2)
    hold on
    plot((1:Nsteps),x_u_series_UKF(1,:),'.-','Color',orange,'LineWidth',2)
    hold on
    plot((1:Nsteps),x_u_series_IPLF(1,:),'.-','Color',black,'LineWidth',2)
    hold on
    plot((1:Nsteps),y,'.-','Color',purple,'LineWidth',2)
    hold on
    plot((1:Nsteps),y_3,'.-','Color',yellow,'LineWidth',2)
    h1 = legend('Truth','Observation','EKF','UKF','IPLF','2nd-order-PLL','3rd-order-PLL','Location','southeast'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',white) %filling colour of legend
    set(h1,'Box','on') %Remove outer frame of legend
    grid on
    box off
    title('KFs vs PLL')
    xlabel('Time step')
    ylabel('Phase[rad]')

    %% Visilization (pre)
x = 1:Nsteps;
M = [x_truth(1,:); y_measure_re; x_u_series_EKF(1,:); x_u_series_UKF(1,:);x_u_series_IPLF(1,:);y';y_3'];
color=[green;brown;blue;orange;black;purple;yellow];
figure;
set(gca,'linewidth',0.4);
set(gca,'GridLineStyle','-.');
set(gca,'GridAlpha',0.4);
set(h1,'Color',white)
set(h1,'Box','on')
grid on
hold on;
lineNames = {'Truth','Observation','EKF','UKF','IPLF','2nd-order-PLL','3rd-order-PLL'};
for i = 1:7 
    plot(x, M(i, :), '.-', 'Color',color(i,:),'LineWidth', 2); 
    legend(lineNames(1:i), 'Location', 'best'); 
    pause;
end

%%
% %% Draw Bode Plot
%     figure(3)
%     bode(H2_s);%check the BW
%     grid on
%     figure(4)
%     bode(H3_s);%check the BW
%     grid on
%% Draw RMSE (part)
    figure(4)
    orange = [1, 0.5, 0];% UKF
    black = [0 0 0];% IPLF
    blue = [0.21 0.35 1];% EKF 
    %yellow = [1, 0.84, 0];% 3 PLL
    %yellow = [0.8, 0.7, 0.1];
    green=[0.4660 0.6740 0.1880];
    purple = [0.54, 0.17, 0.89];% 2 PLL
    white = [1 1 1];
    
    lightgrey = [0.94 0.94 0.94]; % color definition

    plot((1:Nsteps),RMSE_EKF,'--','Color',blue,'LineWidth',2)
    hold on     
    plot((1:Nsteps),RMSE_UKF,'--','Color',orange,'LineWidth',2)
    hold on
    plot((1:Nsteps),RMSE_IPLF,'--','Color',black,'LineWidth',2)


    h1 = legend('rmse(EKF)','rmse(UKF)','rmse(IPLF)','Location','northeast'); 
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
   axis([ 0 Nsteps 0 max(RMSE_EKF)+0.05]) 

% %% Draw prediction (report)
%     orange = [1 0.34 0.20]; %UKF
%     black = [0 0 0]; %IPLF
%     blue = [0.21 0.35 1]; %EKF
%     green = [0.1 0.8 0.5]; %truth
%     purple = [0.75, 0.1, 0.75]; % 2 PLL
%     yellow = [0.8, 0.7, 0.1];% 3 PLL
%     brown = [0.65, 0.16, 0.16]; %observation
%     lightgrey = [0.94 0.94 0.94]; %bg
% 
% 
%     figure;
%     plot((1:Nsteps),x_truth(1,:),'.-','Color',green,'LineWidth',3)
%     hold on    
%     plot((1:Nsteps),y_measure_re,'--','Color',brown,'LineWidth',2)
%     hold on 
%     plot((1:Nsteps),x_u_series_EKF(1,:),'.-','Color',blue,'LineWidth',2)
%     h1 = legend('Truth','Observation','EKF','Location','southeast'); 
%     set(gca,'linewidth',0.4); % thickness of grid
%     set(gca,'GridLineStyle','-.');% type of grid
%     set(gca,'GridAlpha',0.4); % dark of grid
%     set(h1,'Color',white) %filling colour of legend
%     set(h1,'Box','on') %Remove outer frame of legend
%     grid on
%     box off
%     title('EKF')
%     xlabel('Time step')
%     ylabel('Phase[rad]')
% 
%     figure;
%     plot((1:Nsteps),x_truth(1,:),'.-','Color',green,'LineWidth',3)
%     hold on    
%     plot((1:Nsteps),y_measure_re,'--','Color',brown,'LineWidth',2)
%     hold on 
%     plot((1:Nsteps),x_u_series_UKF(1,:),'.-','Color',orange,'LineWidth',2)
%     h1 = legend('Truth','Observation','UKF','Location','southeast'); 
%     set(gca,'linewidth',0.4); % thickness of grid
%     set(gca,'GridLineStyle','-.');% type of grid
%     set(gca,'GridAlpha',0.4); % dark of grid
%     set(h1,'Color',white) %filling colour of legend
%     set(h1,'Box','on') %Remove outer frame of legend
%     grid on
%     box off
%     title('UKF')
%     xlabel('Time step')
%     ylabel('Phase[rad]')
% 
%     figure;
%     plot((1:Nsteps),x_truth(1,:),'.-','Color',green,'LineWidth',3)
%     hold on    
%     plot((1:Nsteps),y_measure_re,'--','Color',brown,'LineWidth',2)
%     hold on 
%     plot((1:Nsteps),x_u_series_IPLF(1,:),'.-','Color',black,'LineWidth',2)
%     h1 = legend('Truth','Observation','IPLF','Location','southeast'); 
%     set(gca,'linewidth',0.4); % thickness of grid
%     set(gca,'GridLineStyle','-.');% type of grid
%     set(gca,'GridAlpha',0.4); % dark of grid
%     set(h1,'Color',white) %filling colour of legend
%     set(h1,'Box','on') %Remove outer frame of legend
%     grid on
%     box off
%     title('IPLF')
%     xlabel('Time step')
%     ylabel('Phase[rad]')
% 
%     figure;
%     plot((1:Nsteps),x_truth(1,:),'.-','Color',green,'LineWidth',3)
%     hold on    
%     plot((1:Nsteps),y_measure_re,'--','Color',brown,'LineWidth',2)
%     hold on 
%     plot((1:Nsteps),y,'.-','Color',purple,'LineWidth',2)
%     h1 = legend('Truth','Observation','2nd-order-PLL','Location','southeast'); 
%     set(gca,'linewidth',0.4); % thickness of grid
%     set(gca,'GridLineStyle','-.');% type of grid
%     set(gca,'GridAlpha',0.4); % dark of grid
%     set(h1,'Color',white) %filling colour of legend
%     set(h1,'Box','on') %Remove outer frame of legend
%     grid on
%     box off
%     title('2nd-order-PLL')
%     xlabel('Time step')
%     ylabel('Phase[rad]')
% 
%     figure;
%     plot((1:Nsteps),x_truth(1,:),'.-','Color',green,'LineWidth',3)
%     hold on    
%     plot((1:Nsteps),y_measure_re,'--','Color',brown,'LineWidth',2)
%     hold on 
%     plot((1:Nsteps),y_3,'.-','Color',yellow,'LineWidth',2)
%     h1 = legend('Truth','Observation','3rd-order-PLL','Location','southeast'); 
%     set(gca,'linewidth',0.4); % thickness of grid
%     set(gca,'GridLineStyle','-.');% type of grid
%     set(gca,'GridAlpha',0.4); % dark of grid
%     set(h1,'Color',white) %filling colour of legend
%     set(h1,'Box','on') %Remove outer frame of legend
%     grid on
%     box off
%     title('3rd-order-PLL')
%     xlabel('Time step')
%     ylabel('Phase[rad]')


