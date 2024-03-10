clear
rng(12,'twister')
%%
%Number of steps
Nsteps = 100;
Nmc=1000; %Number of Monte Carlo runs
f=2;
BW=20*pi;
x_ini=[pi/2,f,0.5]';
P=diag([(pi^2)/3 1 1]);%(pi^2)/3
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix%Q is fixed
R=(pi/3)^2;
H = [1 0 0];%observation matrix
T=0.01;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

%% MC simulation 

independent_R_noise = randn(1, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol=zeros(Nsteps,Nmc);
RMSE_tol_2=zeros(Nsteps,Nmc);
RMSE_tol_3=zeros(Nsteps,Nmc);
for i=1:Nmc
    x_k=x_ini;
    P_k=P;     
    %Measurements
    [x_truth,y_measure_mc] = generate_truth_measurement(Nsteps,x_ini,P,Q,R,H,F,T);
    %SKF
    [x_u_series,RMSE_tol(:,i)] = SKF(Nsteps,x_k,P_k,R,Q,H,F,y_measure_mc,x_truth);  
    %PLL
    [RMSE_tol_2(:,i),RMSE_tol_3(:,i),y,y_3] = PLL_RMSE(y_measure_mc,f,BW,T,Nsteps,x_truth);
end

%% Output mean RMSE
rmse_error_t=sum(RMSE_tol,2)/Nmc;
rmse_pll_2=sum(RMSE_tol_2,2)/Nmc;
rmse_pll_3=sum(RMSE_tol_3,2)/Nmc;
average_RMSE_SKF=sum(rmse_error_t)/Nsteps;
average_RMSE_PLL_2=sum(rmse_pll_2)/Nsteps;
average_RMSE_PLL_3=sum(rmse_pll_3)/Nsteps;
fprintf('RMSE_SKF=%0.5f\n', average_RMSE_SKF);
fprintf('RMSE_PLL_2=%0.5f\n', average_RMSE_PLL_2);
fprintf('RMSE_PLL_3=%0.5f\n', average_RMSE_PLL_3);


%% Draw Tracking
orange = [1, 0.5, 0];
blue = [0.21 0.35 1]; 
yellow = [1, 0.84, 0];
green = [0.1 0.8 0.5];% 3 PLL
purple = [0.54, 0.17, 0.89]; % 2 PLL
brown = [0.65, 0.16, 0.16];
black = [0 0 0];
lightgrey = [0.94 0.94 0.94];
white = [1,1,1];

figure(1)
plot((1:Nsteps),y,'.-','Color', purple,'linewidth',1.5)
hold on
plot((1:Nsteps),y_3,'.-','Color', orange,'linewidth',1.5)
hold on
plot((1:Nsteps),x_u_series(1,:),'.-','Color', blue,'linewidth',1.5)
hold on
plot((1:Nsteps),x_truth(1,:),'.-', 'Color', green,'linewidth',2)
hold on
plot((1:Nsteps),y_measure_mc(1,:),'--','Color', brown,'linewidth',2)

 h1 = legend('2nd-order-PLL','3rd-order-PLL','SKF','Truth','Measurement','Location','northwest'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',white) %filling colour of legend
    set(h1,'Box','on') %Remove outer frame of legend
    grid on
    box off
    title('SKF vs PLLs')
%axis([ 0 Nsteps 0 max(RMSE_2)+0.05]) 
xlabel('Time step')
ylabel('Phase [rad]')


%% Draw RMSE comp
figure(2)
plot (rmse_pll_2,'Color',purple,'linewidth',2)
hold on
plot(rmse_pll_3,'Color',green,'linewidth',2)
hold on
plot(rmse_error_t,'--','Color',black,'linewidth',2)
hold off
 h1 = legend('2nd-order-PLL','3rd-order PLL','SKF','Location','northeast'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',white) %filling colour of legend
    set(h1,'Box','on') %Remove outer frame of legend
    grid on
    box off

title('RMSE comparison')
axis([ 0 Nsteps 0 max(rmse_pll_2)+0.05]) 
xlabel('Time step')
ylabel('RMS phase error [rad]')
% %% Bode plot
% figure(3)
% bode(H2_s);%check the BW
% grid on
% figure(4)
% bode(H3_s);%check the BW
% grid on



    










