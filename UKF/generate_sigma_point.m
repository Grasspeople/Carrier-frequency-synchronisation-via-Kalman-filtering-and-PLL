function [X,W_m,W_c] = generate_sigma_point(x_u, P_u, N_x)
    % N_x - Number of state dimensions
    % W_0 - Weight of the zeroth sigma point
%     lambda=10^(-6)*3-3;
    %sigma points matrix, 2*N_x+1 sigma points
    X=zeros(N_x,2*N_x+1);
    % Calculate the cholesky
    chol_P=chol(3*P_u);
%     chol_P=chol((N_x+lambda)*P_u);
    P_coeff=chol_P; 
    
    % Assign the zeroth sigma point
    X(:,1)=x_u;
    
    % Generate the remaining sigma points
    for i=2:N_x+1

        X(:,i)=x_u+P_coeff(:,i-1);
        X(:,i+N_x)=x_u-P_coeff(:,i-1);
    end
% %***********************************************************************
%   apl=10^(-3);
%   beta=2;
%   W_c = [lambda/(lambda+3)+(1-apl^2+beta); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda))];
%   W_m = [lambda/(lambda+3); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda)); 1/(2*(3+lambda))];
  W_c = [1/3; 1/9; 1/9; 1/9; 1/9; 1/9; 1/9;];
  W_m = [1/3; 1/9; 1/9; 1/9; 1/9; 1/9; 1/9;];%change the w0 to 1/3

end
