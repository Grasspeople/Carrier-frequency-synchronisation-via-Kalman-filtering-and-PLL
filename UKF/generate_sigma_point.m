function [X,W_m,W_c] = generate_sigma_point(x_u, P_u, N_x,lambda)
    % N_x - Number of state dimensions
    % W_0 - Weight of the zeroth sigma point

    %sigma points matrix, 2*N_x+1 sigma points
    X=zeros(N_x,2*N_x+1);
    % Calculate the cholesky
%      chol_P=chol(P_u,'lower');%FIXME:here Pu not positive definite
     chol_P=chol(3*P_u);
%     sqrt_P=sqrtm(P_u);%FIXME:这里的P_u不能正定
    P_coeff=chol_P; %NOTE:为了解决不能chol的问题，我换了sqrtm(P）去尝试
    
    % Assign the zeroth sigma point
    X(:,1)=x_u;
    
    % Generate the remaining sigma points
    for i=2:N_x+1

        X(:,i)=x_u+P_coeff(:,i-1);
        X(:,i+N_x)=x_u-P_coeff(:,i-1);
    end
%***********************************************************************
%   apl=10^(-3);
%   beta=2;
% W_c = [lambda/(lambda+3)+(1-apl^2+beta); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda))];
% W_m = [lambda/(lambda+3); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda))];
  W_c = [1/3; 1/9; 1/9; 1/9; 1/9; 1/9; 1/9;];
  W_m = [1/3; 1/9; 1/9; 1/9; 1/9; 1/9; 1/9;];%change the w0 to 1/3
  
%***********************************************************************


%*********************************    
%     % Weights vector
%     W_m=zeros(2*N_x+1,1);
%     W_m_0=lambda/(N_x+lambda);
%     W_m(1)=W_m_0;
%     
%     % Assign the remaining measurement weights
%     for i =2:2*N_x+1
%         W_m(i)=1/(2*N_x+lambda);
%     end
%     W_c=zeros(2*N_x+1,1);
%     W_c_0=lambda/(N_x+lambda)+3;
%     W_c(1)=W_c_0;
% 
%     % Assign the covariance weights
%     for i =2:2*N_x+1
%         W_c(i)=1/(2*N_x+lambda);
%     end
%*********************************
end
