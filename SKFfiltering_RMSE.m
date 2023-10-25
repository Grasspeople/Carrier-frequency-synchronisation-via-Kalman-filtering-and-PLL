clear
%%
%Number of steps
Nsteps = 100;

x_ini=[pi/2,20,0]';
x_ini_1 = [0,0,0]';
P=diag([(pi^2)/3 1 1]);%(pi^2)/3
Q=diag([0.1 0.1 0.1]);%covariance matrix%Q is fixed
R=diag([3 3 3]);
H = [1, 0, 0; 0, 0, 0; 0, 0, 0];%observation matrix%[1,0,0]?
T=0.05;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure] = generate_truth_measurement(Nsteps,x_ini,P,Q,R,H,F,T);
%%
%Kalman filter
%initialisation
x_k=x_ini_1;
P_k=P;
[x_u_series,sum_error2_squared_t] = SKF(Nsteps,x_k,P_k,R,Q,H,F,y_measure,x_truth);
%%
draw_filtered(Nsteps,y_measure,x_truth,x_u_series)     
%%
%MC simulation
Nmc=10000; %Number of Monte Carlo runs
%Position error (squared)
independent_R_noise = randn(3, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;

for i=1:Nmc
        
    %Measurements
    y_measure_mc=zeros(3,Nsteps);
    y_measure_mc(:,1)=H*x_truth(:,1)+N(:,1);  
    for k=2:Nsteps
        %Measurement
        y_measure_mc(:,k)=H*x_truth(:,k)+N(:,k);
    end
    
    [x_u_series,sum_error2_squared_t] = SKF(Nsteps,x_k,P_k,R,Q,H,F,y_measure,x_truth);  
    
end

rmse_error_t=sqrt(sum_error2_squared_t/Nmc);

figure(2)
plot(rmse_error_t)
ylabel('RMS position error (m)')
xlabel('Time step')
grid on

axis([ 0 Nsteps 0 max(rmse_error_t)+0.05]) 



    










