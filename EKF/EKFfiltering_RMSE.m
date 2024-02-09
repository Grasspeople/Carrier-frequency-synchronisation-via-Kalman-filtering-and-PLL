clear
rng(12,'twister')
%%
%Number of steps
Nsteps =100;

x_ini=[pi/2,20,0]';
P=diag([(pi^2)/3 1 1]);
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);


T=0.01;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure] = generate_truth_measurement(Nsteps,x_ini,Q,R,F);

%%
%Extended Kalman filter
%initialisation
P_k=P;
[x_u_series,RMSE,P_u] = EKF(Nsteps,x_ini,P_k,R,Q,F,y_measure,x_truth);
disp(P_u);

%%
draw_filtered(Nsteps,x_u_series,x_truth)    
%%
%MC simulation
Nmc=1000; %Number of Monte Carlo runs

independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol=zeros(Nsteps,Nmc);

%-----------------------ORIGIN------------------------------------
for i=1:Nmc      
    %Measurements
    [x_truth,y_measure_mc] = generate_truth_measurement(Nsteps,x_ini,Q,R,F);
    [x_u_series,RMSE_tol(:,i)] = EKF(Nsteps,x_ini,P_k,R,Q,F,y_measure_mc,x_truth);
    
end

rmse_error_t=sum(RMSE_tol,2)/Nmc;

figure(2)
plot(rmse_error_t)
ylabel('RMS phase error [rad]')
xlabel('Nsteps')
grid on

axis([ 0 Nsteps 0 max(rmse_error_t)+0.05]) 
averageValue_EKF = mean(rmse_error_t(:))


