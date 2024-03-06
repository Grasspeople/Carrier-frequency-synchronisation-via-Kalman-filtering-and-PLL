clear
rng(12,'twister')
%%
%Number of steps
Nsteps = 100;
% f=1.57542e9/(2*pi); % 1575.42MHz (f的单位是rad)
f=50;
BW=20*pi;
x_ini=[pi/2,f,10]';
%x_ini_1 = [0,0,0]';
P=diag([(pi^2)/3 1 1]);%(pi^2)/3
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix%Q is fixed
R=(pi/3)^2;
H = [1 0 0];%observation matrix
T=1/(5*f);%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure] = generate_truth_measurement(Nsteps,x_ini,P,Q,R,H,F,T);
%% Kalman filter
%initialisation
x_k=x_ini;
P_k=P;
% [x_u_series] = SKF(Nsteps,x_k,P_k,R,Q,H,F,y_measure,x_truth);
% %%
% draw_filtered(Nsteps,y_measure,x_truth,x_u_series)     
%% MC simulation (SKF)
Nmc=1000; %Number of Monte Carlo runs

independent_R_noise = randn(1, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol=zeros(Nsteps,Nmc);
for i=1:Nmc
        
    %Measurements
    [x_truth,y_measure_mc] = generate_truth_measurement(Nsteps,x_ini,P,Q,R,H,F,T);
    [x_u_series,RMSE_tol(:,i)] = SKF(Nsteps,x_k,P_k,R,Q,H,F,y_measure_mc,x_truth);  
end
rmse_error_t=sum(RMSE_tol,2)/Nmc;
average_RMSE_SKF=sum(rmse_error_t)/Nsteps
%% PLL
[average_RMSE_2,RMSE_2,average_RMSE_3,RMSE_3,y,y_3,H2_s,H3_s] = PLL_RMSE(y_measure,f,BW);
%% Draw
orange = [1 0.34 0.20];
blue = [0.21 0.35 1]; 
yellow = [0.1, 0.3, 0.3];
purple = [0.8 0.2 0.5];
black = [0 0 0];
lightgrey = [0.94 0.94 0.94];
white = [1,1,1];
figure(1)
plot((1:Nsteps),y,'b--','linewidth',2)
hold on
plot((1:Nsteps),y_3,'g--','linewidth',2)
hold on
plot((1:Nsteps),x_u_series(1,:),'m.-','linewidth',2)
hold on
plot((1:Nsteps),x_truth(1,:),'r.-','linewidth',2)
hold on
plot((1:Nsteps),y_measure(1,:),'k.-','linewidth',2)

 h1 = legend('2nd-order-PLL','3rd-order-PLL','SKF','Truth','Measurement','Location','northwest'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',white) %filling colour of legend
    set(h1,'Box','on') %Remove outer frame of legend
    grid on
    box off
    title('SKF vs PLLs')



 
figure(2)


plot (RMSE_2,'Color',blue,'linewidth',1)
hold on
plot(rmse_error_t,'Color',orange,'linewidth',1)
 h1 = legend('2nd-order-PLL','SKF','Location','northeast'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',lightgrey) %filling colour of legend
    set(h1,'Box','off') %Remove outer frame of legend
    grid on
    box off

%title('RMSE')
axis([ 0 Nsteps 0 max(RMSE_2)+0.05]) 
xlabel('Time step')
ylabel('RMS phase error [rad]')
%% Bode plot
figure(3)
bode(H2_s);%check the BW
grid on
figure(4)
bode(H3_s);%check the BW
grid on



    










