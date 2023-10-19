clear
orange = [0.66 0.34 0.20]; % 橘色
lightgrey = [0.94 0.94 0.94]; % 淡灰色 % 定义颜色
blue = [0.21 0.35 0.57]; % 蓝色

%%
%Number of steps
Nsteps = 30;

x_ini=[pi/2,0,0]';
P_ini=diag([(pi^2)/3 0.1 0.01]);%(pi^2)/3
Q_ini=diag([0.05 0.01 0.05]);%协方差矩阵
R_ini=diag([0.05 0.01 0.05]);
H = [1, 0, 0; 0, 0, 0; 0, 0, 0];%观测矩阵
chol_P=chol(P_ini)';
%chol_Q=chol(Q_ini)';
chol_R=chol(R_ini)';
independent_P_noise = randn(3, Nsteps);
%independent_Q_noise = randn(3, Nsteps);
independent_R_noise = randn(3, Nsteps);
V = zeros(3,Nsteps);%chol_Q * independent_Q_noise;%过程噪音
N = chol_R * independent_R_noise;%测量噪音

%真实phase/freq/freq dot
x_truth=zeros(3,Nsteps); %theta f fdot
x_truth(:,1) = x_ini + V(:,1);
%测量（只看到phase）
y_measure=zeros(3,Nsteps);
y_measure(:,1)=H*x_truth(:,1)+N(:,1);
%真实相位图
B = (1:Nsteps);
A = vertcat(B, x_truth);
C = vertcat(B, y_measure);

%%
T=0.05;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];
for k=2:Nsteps
    
    %Truth
    x_truth(:,k)=F*x_truth(:,k-1)+V(:,k);
    
    %Measurement
    y_measure(:,k)=H*x_truth(:,k)+N(:,k);
end
%pause

%%
figure(1)
clf
grid on
xlabel('steps')
ylabel('phase [rad]')
plot(A(1,:),x_truth(1,:),'.-b')
hold on
plot(C(1,:),y_measure(1,:),'.-r')
h1 = legend('Truth','Measurements','Location','northwest'); %legend的命令
set(gca,'linewidth',0.4); % 设定网格的粗细
set(gca,'GridLineStyle','-.');% 设定网格的线种类
set(gca,'GridAlpha',0.4); % 设定网格的深浅
set(h1,'Color',lightgrey) %设置legend的填充颜色
set(h1,'Box','off') %去掉legend的外框
box off

%%
%%Kalman filter
%初始化
x_k=x_ini;
P_k=P_ini;
x_u_series=zeros(3,Nsteps);
sum_error2_squared_t=zeros(Nsteps,1);

% figure(2)
% hold on
% plot(x_k(1),x_k(3),'oblack')
% pellipse(x_k([1,3]),P_k([1,3],[1,3]),100,'black',4) %plot 2 sigma region
% grid on
% xlabel('x axis (m)')
% ylabel('y axis (m)')
% axis equal
% axis([-5 40 -10 15])
% hold off


for i=1:Nmc
for k=1:Nsteps
    %Update
    K_gain=P_k*H'/(H*P_k*H'+R_ini);
    x_u=x_k+K_gain*(y_measure(:,k)-H*x_k);
    P_u=P_k-K_gain*H*P_k;
    %Prediction
    x_k=F*x_u;%x_u是前一次的x
    P_k=F*P_u*F';%去掉Q，此模型建模很准确
    
    x_u_series(:,k)=x_u;
    %error加和
    sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u(1))^2;


end    
    figure(2)
%     clf
%     hold on
%     plot(x_u(1),x_u(3),'oblack')
%     pellipse(x_u([1,3]),P_u([1,3],[1,3]),100,'black',4);
%     grid on
%     xlabel('x axis (m)')
%     ylabel('y axis (m)')
%     axis equal
%     axis([-5 40 -10 15])



%%
    xlabel('steps')
    ylabel('phase [rad]')
    
    plot(C(1,:),y_measure(1,:),'.-r')
    hold on
    plot(A(1,:),x_truth(1,:),'.-b')
    hold on
    plot(C(1,:),x_u_series(1,:),'.-g')

    h1 = legend('Truth','Measurements','KF','Location','northwest'); %legend的命令
    set(gca,'linewidth',0.4); % 设定网格的粗细
    set(gca,'GridLineStyle','-.');% 设定网格的线种类
    set(gca,'GridAlpha',0.4); % 设定网格的深浅
    set(h1,'Color',lightgrey) %设置legend的填充颜色
    set(h1,'Box','off') %去掉legend的外框
    box off
    title('KF (green), Measurements (red), True (blue)')

     
%%
%Position error (squared)







