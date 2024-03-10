function [x_u_series,RMSE,cov_pos_j] = IPLF(Nsteps,x_0,P_0,R,Q,F,N_x,x_truth,N_it,y_measure)

%%
x_u_series=zeros(3,Nsteps);
RMSE=zeros(Nsteps,1);
meank=x_0;
Pk=P_0;
%N_it_t=N_it*zeros(1,Nsteps);
for k=1:Nsteps
    
%%   IPLF Update
    mean_pos_j=meank;
    cov_pos_j=Pk;
    for p=1:N_it 
        [X,W] = generate_sigma_point(mean_pos_j, cov_pos_j, N_x);
 
        Y_u_series=zeros(2,2*N_x+1);
        %Value of Y
        for i=1:2*N_x+1
            Y_u = [cos(X(1,i));sin(X(1,i))];
            Y_u_series(:,i)=Y_u;            
        end

        miu=Y_u_series*W;
        S=zeros(2,2);
        C=zeros(3,2);
        for i=1:2*N_x+1
            sub_z_j=Y_u_series(:,i)-miu;
            S=S+W(i)*(sub_z_j*sub_z_j');%predicted covariance of the measurement
            C=C+W(i)*(X(:,i)-mean_pos_j)*sub_z_j';%cross-covariance of the state and the measurement
        end

            %Statistical linear regression parameters
            A_l=C'/cov_pos_j;
            b_l=miu-A_l*mean_pos_j;
            Omega_l=S-A_l*cov_pos_j*A_l';
                    
            %KF (10.15)            miu_j=A_l*x_ini+b_l;
            miu_j=A_l*meank+b_l;
            Psi_j=Pk*A_l';
            S_j=A_l*Pk*A_l'+Omega_l+R;
            
            sub_x=y_measure(:,k)-miu_j;

            mean_update=meank+Psi_j/S_j*sub_x;
            
            var_update=Pk-Psi_j/S_j*Psi_j';
            mean_pos_j=mean_update;
            cov_pos_j=var_update;
    end
        x_u_series(:,k)=mean_update;
         %We sum all errors
        sum_error2_squared_t=zeros(Nsteps);
        sum_error2_squared_t(k)=sum_error2_squared_t(k)+(x_truth(1,k)-x_u_series(1,k))^2;
        RMSE(k)=sqrt(sum_error2_squared_t(k)/k); 
    %% IPLF Prediction
        [X,W] = generate_sigma_point(mean_update, var_update, N_x);% X and W are matrix
        X_hat=F*X;
        var_pred=zeros(3,3);
        x_pred=[0;0;0];
        for j=1:2*N_x+1
           x_pred=x_pred+W(j)*X_hat(:,j);
        end
        for j=1:2*N_x+1
           var_pred=var_pred+W(j)*(X_hat(:,j)-x_pred)*(X_hat(:,j)-x_pred)';
        end
        var_pred=var_pred+Q;        
        meank=x_pred;
        Pk=var_pred;
        Pk=(Pk+Pk')/2;
end  
end

