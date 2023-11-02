clear
%%
%Number of steps
Nsteps = 100;

x_ini=[pi/2,20,0]';
%x_ini_1 = [0,0,0]';
P=diag([(pi^2)/3 1 1]);%(pi^2)/3
Q=diag([0.1 0.1 0.1]);%covariance matrix%Q is fixed
R=(pi/3)^2;
H = [1 0 0];%observation matrix
T=0.05;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure] = generate_truth_measurement(Nsteps,x_ini,P,Q,R,H,F,T);
%%
%Kalman filter
%initialisation
x_k=x_ini;
P_k=P;
[x_u_series] = SKF(Nsteps,x_k,P_k,R,Q,H,F,y_measure,x_truth);
%%
draw_filtered(Nsteps,y_measure,x_truth,x_u_series)     
%%
%MC simulation
Nmc=10000; %Number of Monte Carlo runs

independent_R_noise = randn(1, Nsteps);
chol_R=chol(R)';
N = chol_R * independent_R_noise;
RMSE_tol=zeros(Nsteps,Nmc);
for i=1:Nmc
        
    %Measurements
    [x_truth,y_measure_mc] = generate_truth_measurement(Nsteps,x_ini,P,Q,R,H,F,T);
    [x_u_series,RMSE_tol(:,i)] = SKF(Nsteps,x_k,P_k,R,Q,H,F,y_measure_mc,x_truth);  
    
end

rmse_error_t=sum(RMSE_tol,2)/Nmc;

figure(2)
plot(rmse_error_t)
ylabel('RMS phase error [rad]')
xlabel('Nsteps')
grid on

axis([ 0 Nsteps 0 max(rmse_error_t)+0.05]) 



    










