clear
%%
%Number of steps
Nsteps = 30;

x_ini=[pi/2,0,0]';
P_ini=diag([(pi^2)/3 0.1 0.01]);%(pi^2)/3
Q_ini=diag([0 0 0]);%covariance matrix
R_ini=diag([0.05 0.01 0.05]);
H = [1, 0, 0; 0, 0, 0; 0, 0, 0];%observation matrix
T=0.05;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure] = generate_truth_measurement(Nsteps,x_ini,P_ini,Q_ini,R_ini,H,F,T);
%%
%%Kalman filter
%initialisation
x_k=x_ini;
P_k=P_ini;
[x_u_series] = SKF(Nsteps,x_k,P_k,R_ini,Q_ini,H,F,y_measure);
%%
draw_filtered(Nsteps,y_measure,x_truth,x_u_series)     
%%
%Position error (squared)







