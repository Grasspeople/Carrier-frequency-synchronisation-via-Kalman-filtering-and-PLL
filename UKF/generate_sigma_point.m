function [X,W_m,W_c] = generate_sigma_point(x_ini, P_u, N_x,lambda)
    % N_x - Number of state dimensions
    % W_0 - Weight of the zeroth sigma point

    %sigma points matrix, 2*N_x+1 sigma points
    X=zeros(N_x,2*N_x+1);
    W_m_0=lambda/(N_x+lambda);
    % Calculate the cholesky
    
    chol_P=(chol(P_u))';
    P_coeff=sqrt(N_x/(1-W_m_0))*chol_P;

    
    % Assign the zeroth sigma point
    X(:,1)=x_ini;
    
    % Generate the remaining sigma points
    for i=2:N_x+1

        X(:,i)=x_ini+P_coeff(:,i-1);
        X(:,i+N_x)=x_ini-P_coeff(:,i-1);
    end
    
    % Weights vector
    W_m=zeros(2*N_x+1,1);
    
    W_m(1)=W_m_0;
    
    % Assign the remaining measurement weights
    for i =2:2*N_x+1
        W_m(i)=1/(2*N_x+lambda);
    end

    W_c=zeros(2*N_x+1,1);
    W_c_0=lambda/(N_x+lambda)+2;
    W_c(1)=W_c_0;
    % Assign the covariance weights
    for i =2:2*N_x+1
        W_c(i)=1/(2*N_x+lambda);
    end

end
