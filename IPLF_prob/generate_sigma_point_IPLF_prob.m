function [X,W_m,W_c] = generate_sigma_point_IPLF_prob(x_u, P_u, N_x,lambda)
    % N_x - Number of state dimensions
    W_0=1/3; %- Weight of the zeroth sigma point

    %sigma points matrix, 2*N_x+1 sigma points
    X=zeros(N_x,2*N_x+1);
    % Calculate the cholesky
    chol_P=sqrt(4.5)*chol(P_u);
    P_coeff=chol_P;
    
    % Assign the zeroth sigma point
    X(:,1)=x_u;
    
    % Generate the remaining sigma points
    for i=2:N_x+1

        X(:,i)=x_u+P_coeff(:,i-1);
        X(:,i+N_x)=x_u-P_coeff(:,i-1);
    end
  W_c = [1/3; 1/9; 1/9; 1/9; 1/9; 1/9; 1/9];
  W_m = [1/3; 1/9; 1/9; 1/9; 1/9; 1/9; 1/9];

end
