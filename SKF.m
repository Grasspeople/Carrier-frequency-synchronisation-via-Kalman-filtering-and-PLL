function [x_u_series,sum_error2_squared_t] = SKF(Nsteps,x_k,P_k,R,Q,H,F,y_measure,x_truth)
x_u_series=zeros(3,Nsteps);
sum_error2_squared_t=zeros(Nsteps,1);
for k=1:Nsteps
    %Update
    K_gain=P_k*H'/(H*P_k*H'+R);
    x_u=x_k+K_gain*(y_measure(:,k)-H*x_k);
    P_u=P_k-K_gain*H*P_k;

    %We sum all errors
    sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u(1))^2;
    %Prediction
    x_k=F*x_u;%x_u is former x
    P_k=F*P_u*F'+Q;%remove Q，this model is accurate
    
    x_u_series(:,k)=x_u;
    %RMSE
    %We sum all error
    
end 
end