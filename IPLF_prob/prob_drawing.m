clear
rng(12,'twister')
%% Generate truth and coeff
%Number of steps
Nsteps =100;
x_0=[pi/2,20,0]';
% P_0=diag([42 48 85]);
% Q=diag([56 78 25]);
P_0=diag([(pi^2)/3 1 1]);
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);
N_x=3;
alp=10^(-3);
KAPPA=0; %3-n
lambda=alp^2*(3+KAPPA)-3;% eq.8.74:n=3


T=0.01;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure,Pk] = generate_truth_IPLF_prob(Nsteps,x_0,Q,R,F,P_0);


%% IPLF
% N_it=1;
% Nsteps=1;
% [x_u_series_1,RMSE,Pk_1,var_update] = IPLF_pro(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
% % disp(Pk_1);
% sigma = sqrt(Pk_1(1,1)); % 第一个分量的标准差，s11 是协方差矩阵的第一个对角元素
% mu=x_u_series_1(1,Nsteps);
% x = linspace(mu-3*sigma, mu+3*sigma, 1000);
% pdf_values = normpdf(x, mu, sigma);
% plot(x, pdf_values,'r', 'linewidth', 2);
% hold on
% 
N_it=1;
Nsteps=9;
[x_u_series_3,RMSE,Pk_3,var_update] = IPLF_pro(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
% disp(Pk_2);
sigma = sqrt(Pk_3(1,1));
mu=x_u_series_3(1,Nsteps);
x = linspace(mu-3*sigma, mu+3*sigma, 1000);
pdf_values = normpdf(x, mu, sigma);
plot(x, pdf_values,'b', 'linewidth', 2);
hold on

N_it=5;
Nsteps=9;
[x_u_series_2,RMSE,Pk_2,var_update] = IPLF_pro(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
% disp(Pk_2);
sigma = sqrt(Pk_2(1,1));
mu=x_u_series_2(1,Nsteps);
x = linspace(mu-3*sigma, mu+3*sigma, 1000);
pdf_values = normpdf(x, mu, sigma);
plot(x, pdf_values,'y', 'linewidth', 1);
hold on

N_it=5;
Nsteps=10;
[x_u_series_5,RMSE,Pk_5,var_update] = IPLF_pro(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure);
% disp(Pk_5);
sigma = sqrt(var_update(1,1));
mu=x_u_series_5(1,Nsteps);
x = linspace(mu-3*sigma, mu+3*sigma, 1000);
pdf_values = normpdf(x, mu, sigma);
plot(x, pdf_values,'g');
hold off
xlabel('x1');
ylabel('Probability Density');
% legend('Nsteps=1','Nsteps=10','Nsteps=100')
legend('N_it=1','N_it=5','posterior')




