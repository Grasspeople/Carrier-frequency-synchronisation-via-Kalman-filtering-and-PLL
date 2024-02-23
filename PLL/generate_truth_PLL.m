function [x_truth_phase,y_measure] = generate_truth_PLL(Nsteps,x_ini,Q,R,F)
if Q == zeros(3)
    chol_Q = zeros(3);
else
    chol_Q = chol(Q)';
end
chol_R=chol(R)';
independent_Q_noise = randn(3, Nsteps);
independent_R_noise = randn(2, Nsteps);
V = chol_Q * independent_Q_noise;%chol_Q * independent_Q_noise;%procedure noise
N = chol_R * independent_R_noise;%measurement noise
%%
%true phase/freq/freq dot
x_truth=zeros(3,Nsteps); %theta f fdot
x_truth(:,1) = x_ini + V(:,1);

%measurement（only phase）
y_measure=zeros(2,Nsteps);
h=zeros(2,Nsteps);

h(1,1)=cos(x_truth(1,1));
h(2,1)=sin(x_truth(1,1));
y_measure(:,1)=h(:,1)+N(:,1);
%%
for k=2:Nsteps   
    %Truth
    x_truth(:,k)=F*x_truth(:,k-1)+V(:,k);
    
    %Measurement
    h(1,k)=cos(x_truth(1,k));
    h(2,k)=sin(x_truth(1,k));
    y_measure(:,k)=h(:,k)+N(:,k);
end
x_truth_phase=x_truth(1,:);
end