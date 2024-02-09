function [x_u_series,RMSE,P_u] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure,N_x,x_truth)
x_u_series=zeros(3,Nsteps);
RMSE=zeros(Nsteps,1);
sum_error2_squared_t=zeros(Nsteps);
% x_ini=[0;20;0];
for k=1:Nsteps
    
%%   
    %Update
    [X,W_m,W_c] = generate_sigma_point(x_ini, P, N_x);
    %disp(X)
    
    Y_u_series=zeros(2,2*N_x+1);
    %Value of Y
    for i=1:2*N_x+1
        Y_u = [cos(X(1,i));sin(X(1,i))];
        Y_u_series(:,i)=Y_u;
    end

    miu=zeros(2,1);
    %Predicted mean
    for i=1:2*N_x+1
        miu=miu+W_m(i)*Y_u_series(:,i);  
    end

    S=zeros(2,1);
    %predicted covariance of the measurement
    for i=1:2*N_x+1
        S=S+W_c(i)*(Y_u_series(:,i)-miu)*(Y_u_series(:,i)-miu)';
    end
     S=S+R;

    C=zeros(3,1);
    %cross-covariance of the state and the measurement
    for i=1:2*N_x+1
        C=C+W_c(i)*(X(:,i)-x_ini)*(Y_u_series(:,i)-miu)';
    end

    
    K_gain=C/S;
    x_u=x_ini+K_gain*(y_measure(:,k)-miu);%x_u是mean，x_ini是估计mean
    x_u_series(:,k)=x_u;
    P_u=P-K_gain*S*K_gain';
    %We sum all errors
    sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u_series(1,k))^2;
    RMSE(k)=sqrt(sum_error2_squared_t(k)/k);
%%
    %Prediction
    [X,W_m,W_c] = generate_sigma_point(x_u, P_u, N_x);% X and W are matrix
    X_hat=F*X;
    x_ini=[0,0,0]';
    P=zeros(3);
    for i=1:2*N_x+1
    x_ini=x_ini+W_m(i)*X_hat(:,i);
    end
    for i=1:2*N_x+1
    P=P+W_c(i)*(X_hat(:,i)-x_ini)*(X_hat(:,i)-x_ini)';
    end
    P=P+Q;

%     P=zeros(3);
%     for i=1:2*N_x+1
%     P=P+W_c(i)*X_hat(:,i)*X_hat(:,i)';
%     end
%     x_ini=X_hat*W_m;% 3*7 %7*1
%     P=P-x_ini*x_ini'+Q; 

end 

end