function [x_u_series,RMSE,P_u] = EKF(Nsteps,x_ini,P_k,R,Q,F,y_measure,x_truth, phase,magnitude)
x_u_series=zeros(4,Nsteps);
sum_error2_squared_t=zeros(Nsteps,1);
RMSE=zeros(Nsteps,1);

for k=1:Nsteps
  
    %Update
    h = magnitude(k)*[cos(x_ini(1,1)+phase(k));sin(x_ini(1,1)+phase(k))];
    H = [-sin(x_ini(1,1)+phase(k)) 0 0 -sin(x_ini(1,1)+phase(k));cos(x_ini(1,1)+phase(k)) 0 0 cos(x_ini(1,1)+phase(k))];
    v_k=y_measure(:,k)-h;
    S=H*P_k*H'+R;
    K_gain=P_k*H'/S;
    x_u=x_ini+K_gain*v_k;
    P_u=P_k-K_gain*S*K_gain';

    %Prediction
    x_k=F*x_u;%x_u is former x
    P_k=F*P_u*F'+Q;%remove Q，this model is accurate
    x_ini=x_k;
    x_u_series(:,k)=x_k;

    %We sum all errors
    sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u_series(1,k))^2;
    RMSE(k)=sqrt(sum_error2_squared_t(k)/k);

    
end 

end