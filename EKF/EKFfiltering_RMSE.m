clear
%%
%Number of steps
Nsteps =100;

x_ini=[pi/2,20,0]';
P=diag([(pi^2)/3 1 1]);
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);


T=0.05;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure] = generate_truth_measurement(Nsteps,x_ini,Q,R,F);

%%
%Extended Kalman filter
%initialisation
x_k=x_ini;
P_k=P;

[x_u_series,RMSE] = EKF(Nsteps,x_ini,P_k,R,Q,F,y_measure,x_truth);

%%
draw_filtered(Nsteps,x_u_series,x_truth)    
%%
%MC simulation
Nmc=1000; %Number of Monte Carlo runs

independent_R_noise = randn(2, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol=zeros(Nsteps,Nmc);

% %-------------------------------FIND MAX---------------------------------
% max_RMSE = 0; % 初始化最大RMSE值
% max_i = 0; % 初始化最大RMSE对应的索引i
% x_u_series_all = zeros(Nsteps, Nmc); % 预分配内存，这里的size_of_x_u_series需要您根据实际大小替换
% x_truth_all = zeros(Nsteps, Nmc); % 为x_truth预分配内存
% 
% for i = 1:Nmc
%     %ORIGIN
%     [x_truth,y_measure_mc] = generate_truth_measurement(Nsteps,x_ini,Q,R,F);
%     [x_u_series,RMSE_tol(:,i)] = EKF(Nsteps,x_ini,P_k,R,Q,F,y_measure,x_truth);
%     
%     % 存储每次迭代的x_u_series值
%     x_u_series_all(:, i) = x_u_series(1,:);
%     x_truth_all(:, i) = x_truth(1,:);
%     
%     % 在循环中计算并更新最大RMSE和对应的i
%     current_RMSE = RMSE_tol(:,i);
%     if current_RMSE > max_RMSE
%         max_RMSE = current_RMSE;
%         max_i = i;
%     end
% end
% 
% % 获取最大RMSE对应的数据
% x_u_series_max = x_u_series_all(:, max_i);
% x_truth_max = x_truth_all(:, max_i);
% 
% % 生成图像
% figure(2);
% orange = [1 0.34 0.20];blue = [0.21 0.35 1]; 
% plot(x_u_series_max,'.-','Color',orange);
% hold on
% plot(x_truth_max,'.-','Color',blue);
% 
% xlabel('Time Step');
% ylabel('State Estimate');
% title(sprintf('EKF Performance for Max RMSE at i=%d', max_i));
% legend show;
% hold off;

%-----------------------ORIGIN------------------------------------
rng(0, 'twister')
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


