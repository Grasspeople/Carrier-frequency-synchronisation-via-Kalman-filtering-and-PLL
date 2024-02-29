clear
rng(12,'twister')
%%
%Number of steps
Nsteps =100;
Nmc = 1000;
RMSE_tot=zeros(Nsteps,Nmc);
for j = 1 : Nmc
    x_ini=[pi/2,2000,0]';
    P=diag([(pi^2)/3 1 1]);
    Q=0.1*diag([0.1 0.1 0.1]);%covariance matrix%改大一点
    R=diag([(pi/3)^2 (pi/3)^2]);
    T=0.01;%sampling period
    F=[1 T (T^2)/2; 0 1 T; 0 0 1];
    [x_truth,y_measure] = generate_truth_PLL(Nsteps,x_ini,Q,R,F);
    y=zeros(Nsteps,1);
    y(1)=pi/2;
    x=x_truth;
    xdbuf = pi/2*ones(3,1);
    ydbuf = zeros(2,1);
    ydbuf(1)= pi/2;
    ydbuf(2)= pi/2;
    s=tf('s');
    a=sqrt(2);
    wn=20*pi;%20Hz--40*pi///10Hz--20*pi///5Hz---10*pi
    H2_s=(a*wn*s+wn^2)/(s^2+a*wn*s+wn^2);%open-loop
    H2_z=c2d(H2_s, T, 'foh');
    [num,dem]=tfdata(H2_z,'v');
    sum_error2_squared_t=zeros(Nsteps,1);
    for i = 1:Nsteps %2：Nsteps
        xdbuf = [x(i);xdbuf(1:2)];
        y(i) = (num(1) * xdbuf(1) + num(2) * xdbuf(2) + num(3) * xdbuf(3) - dem(2) * ydbuf(1) - dem(3) * ydbuf(2)) / dem(1);
        ydbuf = [y(i);ydbuf(1)];
    %     %We sum all errors
        sum_error2_squared_t(i)=sum_error2_squared_t(i)+(y(i)-x(i))^2;
        RMSE_tot(i,j)=sqrt(sum_error2_squared_t(i)/i);
    end
end
RMSE=sum(RMSE_tot,2)/Nmc;
average_RMSE=sum(RMSE)/Nsteps

orange = [1 0.34 0.20];
blue = [0.21 0.35 1]; 
lightgrey = [0.94 0.94 0.94];
figure(1)
plot((1:Nsteps),y,'.-','Color',orange)
hold on
plot((1:Nsteps),x,'.-','Color',blue)

 h1 = legend('Input','Output','Location','northwest'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',lightgrey) %filling colour of legend
    set(h1,'Box','off') %Remove outer frame of legend
    grid on
    box off
    title('PLL')

figure(2)
plot(RMSE)
ylabel('RMS phase error [rad]')
xlabel('Nsteps')
grid on
axis([ 0 Nsteps 0 max(RMSE)+0.01]) 
figure(3)
bode(H2_z);%check the BW
grid on



