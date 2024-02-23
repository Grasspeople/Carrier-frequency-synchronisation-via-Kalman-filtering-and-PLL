function [X,W] = sigma_point_IPLF_comp(x_u, P_u, N_x,lambda)
    % N_x - Number of state dimensions
    W0=0.8; %- Weight of the zeroth sigma point

    %sigma points matrix, 2*N_x+1 sigma points
    X=zeros(N_x,2*N_x+1);
    % Calculate the cholesky
    chol_P=sqrt(N_x/(1-W0))*chol(P_u);

    P_coeff=chol_P;
    
    % Assign the zeroth sigma point
    X(:,1)=x_u;
    
    % Generate the remaining sigma points
    for i=2:N_x+1

        X(:,i)=x_u+P_coeff(:,i-1);
        X(:,i+N_x)=x_u-P_coeff(:,i-1);
    end
%   W = [1/3; 1/9; 1/9; 1/9; 1/9; 1/9; 1/9];
  W = [W0; (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x); (1-W0)/(2*N_x)];   %(1-W0)/2Nx=(1-1/3)/2*3=(2/3)/6=1/9

end
