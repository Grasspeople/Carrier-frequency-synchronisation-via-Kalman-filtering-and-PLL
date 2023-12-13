%Author Angel Garcia-Fernandez
clear
randn('seed',12)
rand('seed',12)

mean_ini=[0,2,0,2]';
P_ini=diag([1 0.1 1 0.1]);
chol_ini=chol(P_ini)';

%Nearly constant velocity model
T=0.5;
F=[1 T 0 0;0 1 0 0; 0 0 1 T;0 0 0 1];
sigmaU=0.5;
Q=(sigmaU)^2*[T^3/3 T^2/2 0 0; T^2/2 T 0 0;0 0 T^3/3 T^2/2; 0 0 T^2/2 T];
chol_Q=chol(Q)';


%We measure position
H=[1, 0,0,0;
    0, 0,1,0];
R=0.5*diag([1,1]);
chol_R=chol(R)';

%Number of steps
Nsteps=60;


%Generation of the process
X_truth=zeros(4,Nsteps); %Components 1 and 3 denote position, 2 and 4 velocity
X_truth(:,1)=mean_ini+chol_ini*randn(4,1);
%Measurements
z_series=zeros(2,Nsteps);
z_series(:,1)=H*X_truth(:,1)+chol_R*randn(2,1);


for k=2:Nsteps
    
    %Truth
    X_truth(:,k)=F*X_truth(:,k-1)+chol_Q*randn(4,1);
    plot(X_truth(1,1:k),X_truth(3,1:k),'b')
    
    %Measurement
    z_series(:,k)=H*X_truth(:,k)+chol_R*randn(2,1);
    plot(z_series(1,k),z_series(2,k),'xr')
    %pause(0.5)
end
%pause
plot(z_series(1,:),z_series(2,:),'r')

legend('Truth','Measurements')





%%
%Kalman filter
x_k=mean_ini;
P_k=P_ini;

x_u_series=zeros(4,Nsteps);

figure(2)
hold on
plot(x_k(1),x_k(3),'oblack')
xlabel('x axis (m)')
ylabel('y axis (m)')
axis equal
axis([-5 40 -10 15])
hold off




for k=1:Nsteps
    %Update
    K_gain=P_k*H'/(H*P_k*H'+R);
    x_u=x_k+K_gain*(z_series(:,k)-H*x_k);
    P_u=P_k-K_gain*H*P_k;
    
    x_u_series(:,k)=x_u;
    
    figure(2)
    clf
    hold on
    plot(x_u(1),x_u(3),'oblack')
   
    grid on
    xlabel('x axis (m)')
    ylabel('y axis (m)')
    axis equal
    axis([-5 40 -10 15])
    
    plot(x_u_series(1,1:k),x_u_series(3,1:k),'black')
    plot(X_truth(1,1:k),X_truth(3,1:k),'b')
    plot(z_series(1,1:k),z_series(2,1:k),'xr')
    plot(z_series(1,1:k),z_series(2,1:k),'r')
    hold off
    
    title('KF (black), Measurements (red), True (blue)')

    pause

    
    
    %Prediction
    
    x_k=F*x_u;
    P_k=F*P_u*F'+Q;
    
    
    figure(2)
    clf
    hold on
    plot(x_k(1),x_k(3),'og')
    pellipse(x_k([1,3]),P_k([1,3],[1,3]),100,'g',4);
    grid on
    xlabel('x axis (m)')
    ylabel('y axis (m)')
    axis equal
    axis([-5 40 -10 15])
    
    plot(x_u_series(1,1:k),x_u_series(3,1:k),'black')
    plot(X_truth(1,1:k),X_truth(3,1:k),'b')
    plot(z_series(1,1:k),z_series(2,1:k),'xr')
    plot(z_series(1,1:k),z_series(2,1:k),'r')
    
    hold off
    title('KF (black), Measurements (red), True (blue)')

    pause
   
    
    
end







