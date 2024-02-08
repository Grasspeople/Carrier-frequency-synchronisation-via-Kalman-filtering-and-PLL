function [x_u_series,RMSE] = IPLF(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,lambda,N_it,y_measure)

%%
x_u_series=zeros(3,Nsteps);
RMSE=zeros(Nsteps,1);
sum_error2_squared_t=zeros(Nsteps);
meank=x_0;
Pk=P_0;
%N_it_t=N_it*zeros(1,Nsteps);
for k=1:Nsteps
        miu=zeros(2,1);
        S=zeros(2,2);
        C=zeros(3,2);  
    Y_u_series=zeros(2,2*N_x+1);
%%   IPLF Update
    mean_pos_j=meank;
    cov_pos_j=Pk;
    for p=1:N_it 
        [X,W_m,W_c] = generate_sigma_point_IPLF(mean_pos_j, cov_pos_j, N_x,lambda);
    
        %Value of Y
        for i=1:2*N_x+1
            Y_u = [cos(X(1,i));sin(X(1,i))];
            Y_u_series(:,i)=Y_u;
            miu=miu+W_m(i)*Y_u_series(:,i);%Predicted mean
            S=S+W_c(i)*(Y_u_series(:,i)-miu)*(Y_u_series(:,i)-miu)';%predicted covariance of the measurement
            C=C+W_c(i)*(X(:,i)-mean_pos_j)*(Y_u_series(:,i)-miu)';%cross-covariance of the state and the measurement
        end
            S=S+R;

            %Statistical linear regression parameters
            A_l=C'/cov_pos_j;
            b_l=miu-A_l*mean_pos_j;
            Omega_l=S-A_l*cov_pos_j*A_l';
                    
            %KF (10.15)            miu_j=A_l*x_ini+b_l;
            miu_j=A_l*meank+b_l;

            S_j=A_l*Pk*A_l'+Omega_l;
            K_j=Pk*A_l'/S_j;
            
            sub_x=y_measure(:,k)-miu_j;

            mean_j=meank+K_j*sub_x;
            P_j=Pk-K_j*S_j*K_j';            
            mean_pos_j=mean_j;
            cov_pos_j=P_j; 
    end
        x_u_series(:,k)=mean_j;
         %We sum all errors
        sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u_series(1,k))^2;
        RMSE(k)=sqrt(sum_error2_squared_t(k)/k); 
    %% IPLF Prediction
        [X,W_m,W_c] = generate_sigma_point_IPLF(mean_j, P_j, N_x,lambda);% X and W are matrix
        X_hat=F*X;
        var_pred=zeros(3,3);
        x_pred=[0;0;0];
        for j=1:7
        var_pred=var_pred+W_c(j)*X_hat(:,j)*X_hat(:,j)';
        x_pred=x_pred+W_m(j)*X_hat(:,j);
        end
        var_pred=var_pred+Q;
        
        
        meank=x_pred;
        Pk=var_pred;
        Pk=(Pk+Pk')/2;
end  
end

