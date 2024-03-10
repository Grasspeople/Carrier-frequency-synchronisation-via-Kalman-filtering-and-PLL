function [y,y_3,RMSE_2,RMSE_3] = PLL_RMSE(BW,f,Nsteps,T,y_measure_re,x_truth)

RMSE_2=zeros(Nsteps,1);
RMSE_3=zeros(Nsteps,1);
    y=zeros(Nsteps,1);
    y(1)=pi/2;
    x=y_measure_re(1,:);
    xdbuf = [pi/2-1*T*f;pi/2-2*T*f;0];
    ydbuf = [pi/2-1*T*f;pi/2-2*T*f];
%     xdbuf = (pi/2)*[1;1;1];
%     ydbuf = (pi/2)*[1;1];
    s=tf('s');

    % parameters
    a=sqrt(2); b=1.1; c=2.4;

%     f2 = @(x) x^4-BW^4+4*x^2*BW^2;
%     wn=fzero(f2,10);%61.056Hz 
    wn=BW/((a^2+1)/(4*a^2));
    H2_s=(a*wn*s+wn^2)/(s^2+a*wn*s+wn^2);%open-loop
    H2_z=c2d(H2_s, T, 'foh');
    [num,dem]=tfdata(H2_z,'v');
    sum_error2_squared_t=zeros(Nsteps,1);
    for i = 1:Nsteps %2：Nsteps
        xdbuf = [x(i);xdbuf(1:2)];
        y(i) = (num(1) * xdbuf(1) + num(2) * xdbuf(2) + num(3) * xdbuf(3) - dem(2) * ydbuf(1) - dem(3) * ydbuf(2)) / dem(1);
        ydbuf = [y(i);ydbuf(1)];
    %We sum all errors
        sum_error2_squared_t(i)=sum_error2_squared_t(i)+(y(i)-x_truth(1,i))^2;
        RMSE_2(i)=sqrt(sum_error2_squared_t(i)/i);
    end



%% 3rd-order PLL
    y_3=zeros(Nsteps,1);
    y_3(1)=pi/2;
    x_3=y_measure_re(1,:);
    xdbuf = [pi/2-1*T*f;pi/2-2*T*f;pi/2-3*T*f;0];
    ydbuf = [pi/2-1*T*f;pi/2-2*T*f;pi/2-3*T*f];
%     xdbuf = (pi/2)*[1;1;1;1];
%     ydbuf = (pi/2)*[1;1;1]; 
    s=tf('s');
%     f3 = @(x) (c^2+2*b)*x^2*BW^4+x^6-BW^6+(b^2-2*c)*x^4*BW^2;
%     wn=fzero(f3,20);
    wn=BW/((b*c^2 + (b^2 - c)) / (4*(b*c - 1)));
    H3_s=(c*wn*s^2+b*wn^2*s+wn^3)/(s^3+c*wn*s^2+b*wn^2*s+wn^3);%open-loop
    H3_z=c2d(H3_s, T, 'foh');
    [num_3,dem_3]=tfdata(H3_z,'v');
    sum_error2_squared_t=zeros(Nsteps,1);
    for i = 1:Nsteps %2：Nsteps
        xdbuf = [x_3(i);xdbuf(1:3)];
        y_3(i) = (num_3(1) * xdbuf(1) + num_3(2) * xdbuf(2) + num_3(3) *xdbuf(3) ...
            +num_3(4) * xdbuf(4) - dem_3(2) * ydbuf(1) - dem_3(3) * ydbuf(2)- dem_3(4) * ydbuf(3)) / dem_3(1);
        ydbuf = [y_3(i);ydbuf(1:2)];
    %We sum all errors
        sum_error2_squared_t(i)=sum_error2_squared_t(i)+(y_3(i)-x_truth(1,i))^2;
        RMSE_3(i)=sqrt(sum_error2_squared_t(i)/i);
    end
end






