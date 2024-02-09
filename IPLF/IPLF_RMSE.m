clear
rng(12,'twister')
%% Generate truth and coeff
%Number of steps
Nsteps =100; N_it=1;
x_0=[pi/2,20,0]';
P_0=diag([(pi^2)/3 1 1]);
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);
N_x=3;
alp=10^(-3);
KAPPA=0; %3-n
lambda=alp^2*(3+KAPPA)-3;% eq.8.74:n=3

T=0.01;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure] = generate_truth_IPLF(Nsteps,x_0,Q,R,F);
%% IPLF
[x_u_series,RMSE,cov_pos_j] = IPLF(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
disp(cov_pos_j);
%%
draw_IPLF(Nsteps,x_u_series,x_truth)    
%%
%MC simulation
Nmc=1000; %Number of Monte Carlo runs

independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol=zeros(Nsteps,Nmc);
for i=1:Nmc 
%Measurements
[x_truth,y_measure_mc] = generate_truth_IPLF(Nsteps,x_0,Q,R,F);
[x_u_series,RMSE_tol(:,i)] = IPLF(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
end
rmse_error_t=sum(RMSE_tol,2)/Nmc;
figure(3)
plot(rmse_error_t)
ylabel('RMS phase error [rad]')
xlabel('Nsteps')
grid on
axis([ 0 Nsteps 0 max(rmse_error_t)+0.05]) 

averageValue_IPLF = mean(rmse_error_t(:));
 fprintf('%0.5f\n', averageValue_IPLF);
