function [average_RMSE,RMSE]=draw_filtered(Nsteps,x_u_series_EKF,x_u_series_UKF,x_u_series_IPLF,x_truth,y_measure_re)
    orange = [1 0.34 0.20]; %EKF
    black = [0 0 0]; %UKF
    blue = [0.21 0.35 1]; %truth
    green = [0.1 0.8 0.5]; %observation
    purple = [0.75, 0.1, 0.75]; %1st PLL
    pink = [0.9, 0.5, 0.6]; %IPLF
    lightgrey = [0.94 0.94 0.94]; %bg
    xlabel('steps')
    ylabel('phase [rad]')
    [y,average_RMSE,RMSE] = PLL_RMSE;
    %%    
    figure(1)
    plot((1:Nsteps),x_truth(1,:),'.-','Color',blue)
    hold on     
    plot((1:Nsteps),x_u_series_EKF(1,:),'.-','Color',orange)
    hold on
    plot((1:Nsteps),x_u_series_UKF(1,:),'.-','Color',black)
    hold on
    plot((1:Nsteps),x_u_series_IPLF(1,:),'.-','Color',pink)
    hold on
    plot((1:Nsteps),y_measure_re,'.-','Color',green)
    hold on
    plot((1:Nsteps),y,'.-','Color',purple)
    h1 = legend('Truth','EKF','UKF','Observation','1st-order-PLL','IPLF','Location','northwest'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',lightgrey) %filling colour of legend
    set(h1,'Box','off') %Remove outer frame of legend
    grid on
    box off
    title('KFs vs PLL')
    xlabel('Nsteps')
    ylabel('phase[rad]')

end