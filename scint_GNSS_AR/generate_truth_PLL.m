function [y_measure_re,y_theta,x_truth_phase,y_measure,x_truth] = generate_truth_PLL(Nsteps,x_ini,Q,R,T)

   chol_Q = chol(Q);

   chol_R = chol(R);

independent_Q_noise = randn(4, Nsteps);
independent_R_noise = randn(2, Nsteps);
V = chol_Q * independent_Q_noise;%chol_Q * independent_Q_noise;%procedure noise
N = chol_R * independent_R_noise;%measurement noise
%%
%true phase/freq/freq dot
x_truth=zeros(4,Nsteps); %theta f fdot
x_truth(:,1) = x_ini + V(:,1);

%measurement（only phase）
load('scintDat.mat');
magnitude_full = abs(zkhist);
magnitude = magnitude_full (1:Nsteps);
phase_full = angle(zkhist);
phase = phase_full (1:Nsteps);
% load('scintDat.mat');
% real_full = real(zkhist); % 计算幅度
% realpart = real_full (1:100);
% imag_full = imag(zkhist); % 计算相位，结果以弧度为单位
% imagpart = imag_full (1:100);
y_measure=zeros(2,Nsteps);
h=zeros(2,Nsteps);

h(1,1)=magnitude(1)*cos(x_truth(1,1)+phase(1));
h(2,1)=magnitude(1)*sin(x_truth(1,1)+phase(1));
y_measure(:,1)=h(:,1)+N(:,1);
F=[1 T (T^2)/2; 0 1 T; 0 0 T];
%%
for k=2:Nsteps   
    %Truth
    x_truth(1:3,k)=F*x_truth(1:3,k-1)+V(1:3,k);
    x_truth(4,k)=phase(k)+V(4,k);
    
    %Measurement
    h(1,k)=magnitude(k)*cos(x_truth(1,k)+phase(k));
    h(2,k)=magnitude(k)*sin(x_truth(1,k)+phase(k));
    y_measure(:,k)=h(:,k)+N(:,k);
end
% PLL use it, other don't
complexNumber_y = complex(y_measure(1,:), y_measure(2,:));
y_theta = angle(complexNumber_y)+pi;
complexNumber_N = complex(N(1,:), N(2,:));%都是2pi内的
N_theta=angle(complexNumber_N)+pi;
diff=y_theta-N_theta;
x_truth_phase=x_truth(1,:);
y_measure_re=x_truth_phase+phase(1:Nsteps)'- diff;
% y_measure_re=x_truth_phase+diff;

end