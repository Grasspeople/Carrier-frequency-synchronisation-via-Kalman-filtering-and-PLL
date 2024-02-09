clear
rng(12,'twister')
%% Generate truth and coeff
%Number of steps
Nsteps =100;
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
N_it=1;
N_steps=99;
[x_u_series_1,RMSE_1] = IPLF(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
[f1, xi1] = ksdensity(x_u_series_1(1,:));
plot(xi1, f1, 'r', 'LineWidth', 0.5);
hold on

% N_it=2;
% [x_u_series_2,RMSE_2] = IPLF(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
% [f2, xi2] = ksdensity(x_u_series_2(1,:));
% plot(xi2, f2, 'b');
% hold on

% N_it=3;
% [x_u_series_3,RMSE_3] = IPLF(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
% [f3, xi3] = ksdensity(x_u_series_3(1,:));
% plot(xi3, f3);
% hold on
% 
N_it=10;
[x_u_series_4,RMSE_4] = IPLF(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
[f4, xi4] = ksdensity(x_u_series_4(1,:));
plot(xi4, f4,'b','LineWidth', 0.5);
hold on

[f, xi] = ksdensity(x_truth(1,:));
plot(xi, f, 'g');
% hold on
% 
% % % 绘制PDF
% 
% N_it=5;
% [x_u_series_5,RMSE_5] = IPLF(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
% [f5, xi5] = ksdensity(x_u_series_5(1,:));
% plot(xi5, f5);
hold off

xlabel('theta');
ylabel('Probability Density');


