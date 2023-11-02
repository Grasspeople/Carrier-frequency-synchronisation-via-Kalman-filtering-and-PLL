%This code provides the Kalman filter code and the calculation of the root
%mean square error (RMSE) using Monte Carlo simulation
%Author Angel Garcia-Fernandez
clear
randn('seed',12)
rand('seed',12)

Nmc=100; %Number of Monte Carlo runs

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
Nsteps=30;

%Position error (squared)
sum_error2_squared_t=zeros(Nsteps,1);

%We generate ground truth according to the model
X_truth=zeros(4,Nsteps); %Components 1 and 3 denote position, 2 and 4 velocity
X_truth(:,1)=mean_ini+chol_ini*randn(4,1);
for k=2:Nsteps
    %Truth
    X_truth(:,k)=F*X_truth(:,k-1)+chol_Q*randn(4,1);
end

for i=1:Nmc
        
    %Measurements
    z_series=zeros(2,Nsteps);
    z_series(:,1)=H*X_truth(:,1)+chol_R*randn(2,1);  
    for k=2:Nsteps  
        %Measurement
        z_series(:,k)=H*X_truth(:,k)+chol_R*randn(2,1);
        %plot(z_series(1,k),z_series(2,k),'xr')
        %pause(0.5)
    end
    
    
    
    %%
    %Kalman filter
    x_k=mean_ini;
    P_k=P_ini;
    
    x_u_series=zeros(4,Nsteps);
    
    
    
    for k=1:Nsteps
        %Update
        K_gain=P_k*H'/(H*P_k*H'+R);
        x_u=x_k+K_gain*(z_series(:,k)-H*x_k);
        P_u=P_k-K_gain*H*P_k;        
        x_u_series(:,k)=x_u;
        
        
        %We sum all errors
        sum_error2_squared_t(k)=sum_error2_squared_t(k)+(X_truth(1,k)-x_u(1))^2+(X_truth(3,k)-x_u(3))^2;
        
        %Prediction
        x_k=F*x_u;
        P_k=F*P_u*F'+Q;
             
        
    end
      
    
end

rmse_error_t=sqrt(sum_error2_squared_t/Nmc);

figure(1)
plot(rmse_error_t)
ylabel('RMS position error (m)')
xlabel('Time step')
grid on

axis([ 0 Nsteps 0 max(rmse_error_t)+1]) 
% rmse_error=sqrt(sum(sum_error2_squared_t)/(Nsteps*Nmc));




