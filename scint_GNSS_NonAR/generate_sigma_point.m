function [X,W_m,W_c] = generate_sigma_point(x_u, P_u, N_x)
    % N_x - Number of state dimensions
    % W_0 - Weight of the zeroth sigma point
    W0=0.8;
    %sigma points matrix, 2*N_x+1 sigma points
    X=zeros(N_x,2*N_x+1);
    % Calculate the cholesky
    chol_P=chol((N_x/(1-W0))*P_u);
    P_coeff=chol_P; 
    
    % Assign the zeroth sigma point
    X(:,1)=x_u;
    
    % Generate the remaining sigma points
    for i=2:N_x+1
        X(:,i)=x_u+P_coeff(:,i-1);
        X(:,i+N_x)=x_u-P_coeff(:,i-1);
    end
 
  W_c = [W0; (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x)];
  W_m = [W0; (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x)];%change the w0 to 1/3

end
