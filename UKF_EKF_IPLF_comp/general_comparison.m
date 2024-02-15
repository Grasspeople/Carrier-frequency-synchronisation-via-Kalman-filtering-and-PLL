clear
rng(12,'twister')
%%
%Number of steps
Nsteps = 100; N_it = 5; lambda=10^(-6)*3-3;
x_ini=[pi/2,20,0]';
P=diag([(pi^2)/3 1 1]);
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);
N_x=3;
T=0.01;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];
Nmc=1000; %Number of Monte Carlo runs
independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
UKF_RMSE_tol=zeros(Nsteps,Nmc);
EKF_RMSE_tol=zeros(Nsteps,Nmc);
IPLF_RMSE_tol=zeros(Nsteps,Nmc);
for i=1:Nmc 
% Measurements
[x_truth,y_measure] = truth_comp(Nsteps,x_ini,Q,R,F);
[x_u_series_UKF,UKF_RMSE_tol(:,i)] = UKF_comp(Nsteps,x_ini,P,R,Q,F,y_measure,N_x,x_truth);
[x_u_series_EKF,EKF_RMSE_tol(:,i)] = EKF_comp(Nsteps,x_ini,P,R,Q,F,y_measure,x_truth);
[x_u_series_IPLF,IPLF_RMSE_tol(:,i)] = IPLF_comp(Nsteps,x_ini,P,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);


end
RMSE_UKF=sum(UKF_RMSE_tol,2)/Nmc;
RMSE_EKF=sum(EKF_RMSE_tol,2)/Nmc;
RMSE_IPLF=sum(IPLF_RMSE_tol,2)/Nmc;
figure(1)
plot(RMSE_UKF,LineWidth=1)
hold on
plot(RMSE_EKF,LineWidth=1)
hold on
plot(RMSE_IPLF,LineWidth=1)
hold off
ylabel('RMS phase error [rad]')
xlabel('Nsteps')
grid on
legend('UKF','EKF','IPLF');
% legend([RMSE_UKF, RMSE_EKF, RMSE_IPLF], {'UKF','EKF','IPLF'});
axis([ 0 Nsteps 0 max(RMSE_EKF)+0.05]) 
averageValue_UKF = mean(RMSE_UKF(:));
fprintf('RMSE_UKF=%0.5f\n', averageValue_UKF);
averageValue_EKF = mean(RMSE_EKF(:));
fprintf('RMSE_EKF=%0.5f\n', averageValue_EKF);
averageValue_IPLF = mean(RMSE_IPLF(:));
fprintf('RMSE_IPLF=%0.5f\n', averageValue_IPLF);
