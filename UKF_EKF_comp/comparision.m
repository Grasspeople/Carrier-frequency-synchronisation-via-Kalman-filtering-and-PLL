clear 
n=3;%randam seed
rng(n,'twister');
%%
%Number of steps\

Nsteps = 100;
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
%MC simulation
Nmc=1000; %Number of Monte Carlo runs
independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol_UKF=zeros(Nsteps,Nmc);
RMSE_tol=zeros(Nsteps,Nmc);
x_ini=[pi/2,20,0]';
P_k=P;
for i=1:Nmc     
    %Measurements
    [x_truth,y_measure_mc] = generate_truth_measurement(Nsteps,x_ini,Q,R,F);
    [x_u_series_UKF,RMSE_tol_UKF(:,i)] = UKF(Nsteps,x_ini_UKF,P,R,Q,F,y_measure_mc,N_x,x_truth);
    [x_u_series,RMSE_tol(:,i)] = EKF(Nsteps,x_ini,P_k,R,Q,F,y_measure_mc,x_truth);
    
end
rmse_error_t_UKF=sum(RMSE_tol_UKF,2)/Nmc;
rmse_error_t=sum(RMSE_tol,2)/Nmc;
%%Average
averageValue_UKF = mean(rmse_error_t_UKF(:))
averageValue_EKF = mean(rmse_error_t(:))





