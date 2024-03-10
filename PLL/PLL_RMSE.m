clear
rng(12,'twister')
%%
%Number of steps
Nsteps =100;
Nmc = 1000;
f=2000;
BW=100*pi;%20Hz--40*pi///10Hz--20*pi///5Hz---10*pi
RMSE_tot=zeros(Nsteps,Nmc);
RMSE_tot_3=zeros(Nsteps,Nmc);
%% 2nd-order PLL
for j = 1 : Nmc
    x_ini=[pi/2,f,0]';
    P=diag([(pi^2)/3 1 1]);
    Q=0.1*diag([0.1 0.1 0.1]);%covariance matrix%改大一点
    R=diag([(pi/3)^2 (pi/3)^2]);
    T=1/(5*f);%sampling period
    F=[1 T (T^2)/2; 0 1 T; 0 0 1];
    [y_measure_re,y_theta,x_truth_phase,y_measure,x_truth] = generate_truth_PLL(Nsteps,x_ini,Q,R,F);
    y=zeros(Nsteps,1);
    y(1)=pi/2;
    x=y_measure_re;
    xdbuf = 0*ones(4,1);
    ydbuf = zeros(3,1);
%     xdbuf = (pi/2)*ones(3,1);
%     ydbuf = zeros(2,1);
%     ydbuf(1)= pi/2;
%     ydbuf(2)= pi/2;
    s=tf('s');
    a=sqrt(2);
    wn=BW/((a^2+1)/(4*a^2));
    H2_s=(a*wn*s+wn^2)/(s^2+a*wn*s+wn^2);%open-loop
    H2_z=c2d(H2_s, T, 'foh');
    [num,dem]=tfdata(H2_z,'v');
    sum_error2_squared_t=zeros(Nsteps,1);
    for i = 1:Nsteps %2：Nsteps
        xdbuf = [x(i);xdbuf(1:2)];
        y(i) = (num(1) * xdbuf(1) + num(2) * xdbuf(2) + num(3) * xdbuf(3) - dem(2) * ydbuf(1) - dem(3) * ydbuf(2)) / dem(1);
        ydbuf = [y(i);ydbuf(1)];
    %We sum all errors
        sum_error2_squared_t(i)=sum_error2_squared_t(i)+(y(i)-x(i))^2;
        RMSE_tot(i,j)=sqrt(sum_error2_squared_t(i)/i);
    end
end
RMSE=sum(RMSE_tot,2)/Nmc;
average_RMSE=sum(RMSE)/Nsteps

b=1.1;c=2.4;
%% 3rd-order PLL
for j = 1 : Nmc
    x_ini=[pi/2,f,0]';
    P=diag([(pi^2)/3 1 1]);
    Q=0.1*diag([0.1 0.1 0.1]);%covariance matrix%改大一点
    R=diag([(pi/3)^2 (pi/3)^2]);
    T=1/(5*f);%sampling period
    F=[1 T (T^2)/2; 0 1 T; 0 0 1];
    [y_measure_re,y_theta,x_truth_phase,y_measure,x_truth] = generate_truth_PLL(Nsteps,x_ini,Q,R,F);
    y_3=zeros(Nsteps,1);
    y_3(1)=pi/2;
    x_3=y_measure_re;
    xdbuf = 0*ones(4,1);
    ydbuf = zeros(3,1);
%     ydbuf(1)= pi/2;
%     xdbuf = [pi/2-T*2*pi*f;pi/2-2*T*2*pi*f;pi/2-3*T*2*pi*f;0];
%     ydbuf = [pi/2-T*2*pi*f;pi/2-2*T*2*pi*f;pi/2-3*T*2*pi*f];
    s=tf('s');
    a=sqrt(2);
    wn=BW/((b*c^2 + (b^2 - c)) / (4*(b*c - 1)));
    H3_s=(c*wn*s^2+b*wn^2*s+wn^3)/(s^3+c*wn*s^2+b*wn^2*s+wn^3);%open-loop
    

    H3_z=c2d(H3_s, T, 'foh');
    [num_3,dem_3]=tfdata(H3_z,'v');
    sum_error2_squared_t=zeros(Nsteps,1);
    for i = 1:Nsteps %2：Nsteps
        xdbuf = [x_3(i);xdbuf(1:3)];
        y_3(i) = (num_3(1) * xdbuf(1) + num_3(2) * xdbuf(2) + num_3(3) *xdbuf(3) ...
            +num_3(4) * xdbuf(4) - dem_3(2) * ydbuf(1) - dem_3(3) * ydbuf(2)- dem_3(4) * ydbuf(3)) / dem_3(1);
        ydbuf = [y_3(i);ydbuf(1:2)];
    %We sum all errors
        sum_error2_squared_t(i)=sum_error2_squared_t(i)+(y_3(i)-x_3(i))^2;
        RMSE_tot_3(i,j)=sqrt(sum_error2_squared_t(i)/i);
    end
end
RMSE_3=sum(RMSE_tot_3,2)/Nmc;
average_RMSE_3=sum(RMSE_3)/Nsteps

%% Draw
orange = [1 0.34 0.20];
blue = [0.21 0.35 1]; 
green = [0.21 0.7 0.5];
purple = [0.8 0.2 0.5];
lightgrey = [0.94 0.94 0.94];
figure(1)
plot((1:Nsteps),y,'.-','Color',orange)
hold on
plot((1:Nsteps),y_3,'.-','Color',purple)
hold on
plot((1:Nsteps),x_truth_phase,'.-','Color',blue)
hold on
plot((1:Nsteps),y_measure_re,'.-','Color',green)

 h1 = legend('2nd-order-PLL','3rd-order-PLL','Truth','Measurement','Location','northwest'); 
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
xlabel('Time step')
grid on
axis([ 0 Nsteps 0 max(RMSE)+0.01]) 
figure(3)
bode(H2_s);%check the BW
grid on
figure(4)
bode(H3_s);%check the BW
grid on



