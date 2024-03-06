function [average_RMSE_2,average_RMSE_3,RMSE_2,RMSE_3]=draw_filtered(Nsteps,x_u_series_EKF,x_u_series_UKF,x_u_series_IPLF,x_truth,y_measure_re,BW,f)
    
    orange = [1, 0.5, 0]; 
    black = [0 0 0];
    blue = [0.21 0.35 1]; 
    green = [0.2, 0.8, 0.2]; 
    purple = [0.54, 0.17, 0.89];
%     brown = [0.5, 0.16, 0.16];
    yellow = [1, 0.84, 0];% 3 PLL
    



    orange = [1 0.34 0.20]; %UKF
    black = [0 0 0]; %IPLF
    blue = [0.21 0.35 1]; %EKF
    green = [0.1 0.8 0.5]; %truth
    purple = [0.75, 0.1, 0.75]; % 2 PLL
    yellow = [0.8, 0.7, 0.1];% 3 PLL
    brown = [0.65, 0.16, 0.16]; %observation
    lightgrey = [0.94 0.94 0.94]; %bg
    xlabel('steps')
    ylabel('phase [rad]')
    [y,y_3,average_RMSE_2,average_RMSE_3,RMSE_2,RMSE_3] = PLL_RMSE(BW,f);
    %%    
    figure(1)
    plot((1:Nsteps),x_truth(1,:),'.-','Color',green,'LineWidth',3)
    hold on     
    plot((1:Nsteps),x_u_series_EKF(1,:),'.-','Color',blue,'LineWidth',1)
    hold on
    plot((1:Nsteps),x_u_series_UKF(1,:),'.-','Color',orange,'LineWidth',1)
    hold on
    plot((1:Nsteps),x_u_series_IPLF(1,:),'.-','Color',black,'LineWidth',1)
    hold on
    plot((1:Nsteps),y_measure_re,'--','Color',brown,'LineWidth',2)
    hold on
    plot((1:Nsteps),y,'.-','Color',purple,'LineWidth',1)
    hold on
    plot((1:Nsteps),y_3,'.-','Color',yellow,'LineWidth',1)
    h1 = legend('Truth','EKF','UKF','IPLF','Observation','2nd-order-PLL','3rd-order-PLL','Location','northwest'); 
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