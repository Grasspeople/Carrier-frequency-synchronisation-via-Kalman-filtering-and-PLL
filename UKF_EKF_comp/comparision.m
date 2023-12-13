clear
rng(12,'twister');
orange = [1 0.34 0.20]; 
red = [0.5 0 1]; % color definition
blue = [0.21 0.35 1]; 
%%
%Number of steps
Nsteps =100;
x_ini_UKF=[pi/2,20,0]';
P=diag([(pi^2)/3 1 1]);
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);
N_x=3;
lambda=10^(-5);
T=0.05;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];
%%
%UKF
[x_truth,y_measure] = generate_truth_measurement(Nsteps,x_ini_UKF,Q,R,F);
[x_u_series_UKF] = UKF(Nsteps,x_ini_UKF,P,R,Q,F,y_measure,N_x,x_truth);
%MC simulation
Nmc=1000; %Number of Monte Carlo runs
independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol=zeros(Nsteps,Nmc);
for i=1:Nmc 
%Measurements
[x_truth,y_measure_mc] = generate_truth_UKF(Nsteps,x_ini_UKF,Q,R,F);
[x_u_series_UKF,RMSE_tol(:,i)] = UKF(Nsteps,x_ini_UKF,P,R,Q,F,y_measure_mc,N_x,x_truth);
end
rmse_error_t_UKF=sum(RMSE_tol,2)/Nmc;
%%Average
averageValue_UKF = mean(rmse_error_t_UKF(:))
%%
%EKF
x_ini=[pi/2,20,0]';
x_k=x_ini;
P_k=P;
[x_u_series] = EKF(Nsteps,x_ini,P_k,R,Q,F,y_measure,x_truth);
%MC simulation
independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
RMSE_tol=zeros(Nsteps,Nmc);
for i=1:Nmc 
%Measurements
[x_truth,y_measure_mc] = generate_truth_UKF(Nsteps,x_ini,Q,R,F);
[x_u_series_UKF,RMSE_tol(:,i)] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure_mc,N_x,x_truth);
end
rmse_error_t=sum(RMSE_tol,2)/Nmc;
%%Average
averageValue_EKF = mean(rmse_error_t(:))

figure(1)
plot((1:Nsteps),x_truth(1,:),'.-','Color',blue)
hold on     
plot((1:Nsteps),x_u_series(1,:),'.-','Color',orange)
hold on
plot((1:Nsteps),x_u_series_UKF(1,:),'.-','Color',red)
title('Tracking Comparision')
legend('truth','EKF','UKF','Location','northwest'); 
xlabel('Nsteps')
ylabel('phase[rad]')





