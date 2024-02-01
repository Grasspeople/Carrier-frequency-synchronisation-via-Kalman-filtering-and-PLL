clear
rng(12,'twister')
%% Generate truth and coeff
%Number of steps
Nsteps =100; N_it=10;

mean_pre=[pi/2,20,0]';
cov_pre=diag([(pi^2)/3 1 1]);
Q=0.01*diag([0.1 0.1 0.1]);%covariance matrix
R=diag([(pi/3)^2 (pi/3)^2]);
N_x=3;
alp=10^(-3);
KAPPA=0; %3-n
lambda=alp^2*(3+KAPPA)-3;% eq.8.74:n=3
disp(lambda)

T=0.05;%sampling period
F=[1 T (T^2)/2; 0 1 T; 0 0 1];

[x_truth,y_measure] = generate_truth_IPLF(Nsteps,mean_pre,Q,R,F);
%%
x_u_series=zeros(3,Nsteps);
x_u_series(:,1)=[pi/2,20,0]';
RMSE=zeros(Nsteps,1);
sum_error2_squared_t=zeros(Nsteps);
N_it_t=N_it*zeros(1,Nsteps);
for k=1:Nsteps
    
    Y_u_series=zeros(2,2*N_x+1);
    miu=zeros(2,1);
    S=zeros(2,1);
    C=zeros(3,1);
%%   IPLF Update
    for p=1:N_it 
        [X,W_m,W_c] = generate_sigma_point_IPLF(mean_pre, cov_pre, N_x,lambda);
        %disp(X)
    
        %Value of Y
        for i=1:2*N_x+1
            Y_u = [cos(X(1,i));sin(X(1,i))];
            Y_u_series(:,i)=Y_u;
        end
    
        %Predicted mean
        for i=1:2*N_x+1
            miu=miu+W_m(i)*Y_u_series(:,i);  
        end

        %predicted covariance of the measurement
        for i=1:2*N_x+1
            S=S+W_c(i)*(Y_u_series(:,i)-miu)*(Y_u_series(:,i)-miu)';
        end
        S=S+R;

        %cross-covariance of the state and the measurement
        for i=1:2*N_x+1
            C=C+W_c(i)*(X(:,i)-mean_pre)*(Y_u_series(:,i)-miu)';
        end


        %计算miu,C,S和UKF一样，以下引入新的变量：A,b,Omega
            %Statistical linear regression parameters
            A_l=C'/cov_pre;
            b_l=miu-A_l*mean_pre;
            Omega_l=S-A_l*cov_pre*A_l';
                    
            %KF (10.15)            miu_j=A_l*x_ini+b_l;
            miu_j=A_l*mean_pre+b_l;
            S_j=A_l*cov_pre*A_l'+Omega_l;
            K_gain_j=cov_pre*A_l'/S_j;

            sub_x=x_truth(1,k)-miu_j;

            mean_j=mean_pre+K_gain_j*sub_x;
            P_j=cov_pre-K_gain_j*S_j*K_gain_j';
            
            N_it_t(k)=p;
            mean_pre=mean_j;%FIXME:m_k只是书上记法，和代码里不一样
            cov_pre=P_j; %FIXME:P_k只是书上记法，和代码里不一样

    end
    x_u_series(:,k)=mean_pre;
%% IPLF Prediction
    [X,W_m,W_c] = generate_sigma_point_IPLF(mean_pre, cov_pre, N_x,lambda);% X and W are matrix
    X_hat=F*X;
    mean_pre=[0,0,0]';%这是m-
    for i=1:2*N_x+1
    mean_pre=mean_pre+W_m(i)*X_hat(:,i);%这是m-，已计算出
    end

    cov_pre=zeros(3);
    for i=1:2*N_x+1
    cov_pre=cov_pre+W_c(i)*(X_hat(:,i)-mean_pre)*(X_hat(:,i)-mean_pre)';
    end
    cov_pre=cov_pre+Q;%此时P计算出
    cov_pre=(cov_pre+cov_pre')/2;
    
    

%     %We sum all errors
%     sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u_series(1,k))^2;
%     RMSE(k)=sqrt(sum_error2_squared_t(k)/k);

    
end 
draw_IPLF(Nsteps,x_u_series,x_truth)