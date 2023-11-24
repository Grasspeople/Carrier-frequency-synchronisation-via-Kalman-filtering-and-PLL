function [X,W] = generate_sigma_point(x_u, P_u, N_x, W_0)
    % N_x - Number of state dimensions
    % W_0 - Weight of the zeroth sigma point

    %sigma points matrix, 2*N_x+1 sigma points
    X=zeros(N_x,2*N_x+1);
    
    % Calculate the cholesky
    P_coeff=N_x/(1-W_0)*P_u;
    chol_P=chol(P_coeff,'lower');
    
    % Assign the zeroth sigma point
    X(:,1)=x_u;
    
    % Generate the remaining sigma points
    for i=2:N_x+1
        X(:,i)=x_u+chol_P(:,i);
        X(:,i+N_x)=x_u-chol_P(:,i);
    end
    
    % Weights vector
    W=zeros(2*N_x+1,1);
    W(1)=W_0;
    
    % Assign the remaining weights
    for i =2:2*N_x+1
        W(i)=(1-W_0)/(2*N_x);
    end
end
