function [x_u_series,RMSE] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure,N_x,x_truth)
x_u_series=zeros(3,Nsteps);
%sum_error2_squared_t=zeros(Nsteps,1);
x_u_series(:,1)=[pi/2,20,0]';
RMSE=zeros(Nsteps,1);
lambda=0;
sum_error2_squared_t=zeros(Nsteps);

for k=1:Nsteps
    
    Y_u_series=zeros(2,2*N_x+1);
    miu=zeros(2,1);
    S=zeros(2,1);
    C=zeros(3,1);
%%   
    %Update
    [X,W_m,W_c] = generate_sigma_point(x_ini, P, N_x,lambda);
    
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
        C=C+W_c(i)*(X(:,i)-x_ini)*(Y_u_series(:,i)-miu)';
    end

    
    K_gain=C/S;
    x_u=x_ini+K_gain*(y_measure(:,k)-miu);
    P_u=P-K_gain*S*K_gain';
    %P_u(abs(P_u) < 1e-15) = 0;%solution
    x_u_series(:,k)=x_u;

%%
    %Prediction
    [X,W_m,W_c] = generate_sigma_point(x_u, P_u, N_x,lambda);% X and W are matrix
    
    X_hat=F*X;
    x_ini=[0;0;0];
    for i=1:2*N_x+1
    x_ini=x_ini+W_m(i)*X_hat(:,i);
    end

    P=zeros(3);
    for i=1:2*N_x+1
    P=P+W_c(i)*(X_hat(:,i)-x_ini)*(X_hat(:,i)-x_ini)';
    end
    P=P+Q;

    

    %We sum all errors
    sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u_series(1,k))^2;
    RMSE(k)=sqrt(sum_error2_squared_t(k)/k);

    
end 

end