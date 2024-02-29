function [y,average_RMSE,RMSE] = PLL_RMSE
rng(12,'twister')
%%
%Number of steps
Nsteps =100;
Nmc = 1000;
load('scintDat.mat');
magnitude_full = abs(zkhist);
magnitude = magnitude_full (1:100);
phase_full = angle(zkhist);
phase = phase_full (1:100);
f=2000; % frequency of carrier signal 
RMSE_tot=zeros(Nsteps,Nmc);
for j = 1 : Nmc
    x_ini=[pi/2, f, 0, phase(1)]';
    P=diag([(pi^2)/3 1 1 0.08]); % 0.08 is from Valis paper
    Q=0.1*diag([0.1 0.1 0.1 0.1]);%covariance matrix
    R=diag([(pi/3)^2 (pi/3)^2]);
    T=1/(5*f);%sampling period
    F=[1 T (T^2)/2; 0 1 T; 0 0 1];
    [y_measure_re,y_theta,x_truth_phase] = generate_truth_PLL(Nsteps,x_ini,Q,R,T);
    y=zeros(Nsteps,1);
    y(1)=pi/2;
    x=y_measure_re;
    xdbuf = pi/2*ones(3,1);
    ydbuf = zeros(2,1);
    ydbuf(1)= pi/2;
    ydbuf(2)= pi/2;
    s=tf('s');
    a=sqrt(2);
    wn=20*pi;%20Hz--40*pi///10Hz--20*pi///5Hz---10*pi (10Hz)
    H2_s=(a*wn*s+wn^2)/(s^2+a*wn*s+wn^2);%open-loop
    H2_z=c2d(H2_s, T, 'foh');
    [num,dem]=tfdata(H2_z,'v');
    sum_error2_squared_t=zeros(Nsteps,1);
    for i = 1:Nsteps %2：Nsteps
        xdbuf = [x(i);xdbuf(1:2)];
        y(i) = (num(1) * xdbuf(1) + num(2) * xdbuf(2) + num(3) * xdbuf(3) - dem(2) * ydbuf(1) - dem(3) * ydbuf(2)) / dem(1);
        ydbuf = [y(i);ydbuf(1)];
    %     %We sum all errors
        sum_error2_squared_t(i)=sum_error2_squared_t(i)+(y(i)-x_truth_phase(i))^2;
        RMSE_tot(i,j)=sqrt(sum_error2_squared_t(i)/i);
    end
end
RMSE=sum(RMSE_tot,2)/Nmc;
average_RMSE=sum(RMSE)/Nsteps
end




