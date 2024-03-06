function [x_u_series,RMSE,P_u] = UKF(Nsteps,x_ini,P,R,Q,F,y_measure,N_x,x_truth,magnitude)
x_u_series=zeros(4,Nsteps);
RMSE=zeros(Nsteps,1);

x_pred=x_ini;
P_pred=P;
ob_pred_point_series=zeros(2,2*N_x+1);
for k=1:Nsteps   
    
     %Update
    [X,W_m,W_c] = generate_sigma_point(x_pred, P_pred, N_x);
    for i=1:2*N_x+1
        ob_pred_point=magnitude(k)*[cos(X(1,i)+X(4,i));sin(X(1,i)+X(4,i))];
        ob_pred_point_series(:,i)=ob_pred_point;
    end
    
    S=zeros(2,2);
    ob_pred = W_m(i) * ob_pred_point_series(:,i);  
    for i=1:2*N_x+1
        S = S + W_c(i) * (ob_pred_point_series(:,i) - ob_pred) * (ob_pred_point_series(:,i) - ob_pred)';
    end
    S=S+R;

    C=zeros(4,2);
    for i=1:2*N_x+1
        C=C+W_c(i) * (X(:,i)-x_pred)*(ob_pred_point_series(:,i)-ob_pred)';
    end 
    K_gain=C/S;
    minus=y_measure(:,k)-ob_pred;
    x_u=x_pred+K_gain*minus;
    x_u_series(:,k)=x_u;
    P_u=P_pred-K_gain*S*K_gain';
    %We sum all errors
    sum_error2_squared_t=zeros(Nsteps);
    sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u_series(1,k))^2;
    RMSE(k)=sqrt(sum_error2_squared_t(k)/k);
    
    %Prediction
    [X,W_m,W_c] = generate_sigma_point(x_u, P_u, N_x);
    X_hat=F*X;
    x_pred=X_hat*W_m;
    P_pred=zeros(4,4);
    for i=1:2 * N_x + 1
    P_pred = P_pred + W_c(i) * (X_hat(:,i) - x_pred) * (X_hat(:,i) - x_pred)';
    end
    P_pred=P_pred + Q;

end 
end