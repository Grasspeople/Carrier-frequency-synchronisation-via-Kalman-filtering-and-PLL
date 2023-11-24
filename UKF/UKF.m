function [x_u_series,RMSE] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure,x_truth,N_x)
x_u_series=zeros(3,Nsteps);
%sum_error2_squared_t=zeros(Nsteps,1);
x_u_series(:,1)=x_ini;
RMSE=zeros(Nsteps,1);


for k=1:Nsteps
    
    Y_u_series=zeros(2,2*N_x+1);
    miu=zeros(2,1);
    S=zeros(2,1);

    %Update
    %Value of Y
    for i=1:2*N_x+1
        Y_u = [cos(x_ini(1));sin(x_ini(1))];
        Y_u_series(:,i)=Y_u;
    end
    
    %Predicted mean
    for i=1:2*N_x+1
        miu=miu+W(i)*Y_u_series(:,i);   %W(i) unrecognized, need W_ini(?)
    end

    %predicted covariance of the measurement
    for i=1:2*N_x+1
        S=S+W(i)*(Y_u_series(:,i)-miu)*(Y_u_series(:,i)-miu)'+R;
    end

    %cross-covariance of the state and the measurement
    for i=1:2*N_x+1
        C=C+W(i)*(X(:,i)-x_ini)*(Y_u_series(:,i)-miu)';
    end

    K_gain=C/S;
    x_u=x_ini+K_gain*(y_measure(:,Nsteps)-miu);
    P_u=P-K_gain*S*K_gain';
    x_u_series(:,k)=x_u;

    [X,W] = generate_sigma_point(x_u, P_u, N_x, W_0);% X and W are matrix

    %Prediction
    X_hat=F*X;
    for i=1:2*N_x+1
    x_ini=x_ini+W(i)*X_hat;
    end


    for i=1:2*N_x+1
    P=P+W(i)*(X_hat-x_ini)*(X_hat-x_ini)'+Q;
    end

    

%     %We sum all errors
%     sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u_series(1,k))^2;
%     RMSE(k)=sum_error2_squared_t(k)/k;

    
end 

end