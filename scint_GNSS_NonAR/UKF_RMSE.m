clear
rng(12,'twister')
%%
%Number of steps
Nsteps =100;

x_ini=[pi/2,2000,0]';
P=diag([(pi^2)/3 1 1]);
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);
N_x=3;
T=0.01;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure] = generate_truth_UKF(Nsteps,x_ini,Q,R,F);

%%
%Unscented Kalman filter
[x_u_series,RMSE,P_u] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure,N_x,x_truth);

%%
draw_UKF(Nsteps,x_u_series,x_truth)    
%%
% MC simulation
Nmc=1000; %Number of Monte Carlo runs
independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol=zeros(Nsteps,Nmc);
for i=1:Nmc 
% Measurements
[x_truth,y_measure_mc] = generate_truth_UKF(Nsteps,x_ini,Q,R,F);
[x_u_series_UKF,RMSE_tol(:,i)] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure_mc,N_x,x_truth);
end
rmse_error_t=sum(RMSE_tol,2)/Nmc;
figure(3)
plot(rmse_error_t)
ylabel('RMS phase error [rad]')
xlabel('Nsteps')
grid on
axis([ 0 Nsteps 0 max(rmse_error_t)+0.05]) 
averageValue_UKF = mean(rmse_error_t(:));
fprintf('%0.5f\n', averageValue_UKF);
