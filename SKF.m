function [x_u_series] = SKF(Nsteps,x_k,P_k,R_ini,Q_ini,H,F,y_measure)
x_u_series=zeros(3,Nsteps);
for k=1:Nsteps
    %Update
    K_gain=P_k*H'/(H*P_k*H'+R_ini);
    x_u=x_k+K_gain*(y_measure(:,k)-H*x_k);
    P_u=P_k-K_gain*H*P_k;
    %Prediction
    x_k=F*x_u;%x_u is former x
    P_k=F*P_u*F'+Q_ini;%remove Q，this model is accurate
    
    x_u_series(:,k)=x_u;


end    
end