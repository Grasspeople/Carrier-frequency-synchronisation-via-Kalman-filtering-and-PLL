function draw_filtered(Nsteps,x_u_series_EKF,x_u_series_UKF,x_u_series_IPLF,x_truth,y_measure_re,y,y_3)
    orange = [1 0.34 0.20]; %UKF
    black = [0 0 0]; %IPLF
    blue = [0.21 0.35 1]; %EKF
    green = [0.1 0.8 0.5]; %truth
    purple = [0.75, 0.1, 0.75]; % 2 PLL
    yellow = [0.8, 0.7, 0.1];% 3 PLL
    white=[1 1 1];
    brown = [0.65, 0.16, 0.16]; %observation
    lightgrey = [0.94 0.94 0.94]; %bg

    %%    
    figure(1)
    plot((1:Nsteps),x_truth(1,:),'.-','Color',green,'LineWidth',3)
    hold on     
    plot((1:Nsteps),x_u_series_EKF(1,:),'.-','Color',blue,'LineWidth',1.5)
    hold on
    plot((1:Nsteps),x_u_series_UKF(1,:),'.-','Color',orange,'LineWidth',1.5)
    hold on
    plot((1:Nsteps),x_u_series_IPLF(1,:),'.-','Color',black,'LineWidth',1.5)
    hold on
    plot((1:Nsteps),y_measure_re,'--','Color',brown,'LineWidth',2)
    hold on
    plot((1:Nsteps),y,'.-','Color',purple,'LineWidth',1.5)
    hold on
    plot((1:Nsteps),y_3,'.-','Color',yellow,'LineWidth',1.5)
    h1 = legend('Truth','EKF','UKF','IPLF','Observation','2nd-order-PLL','3rd-order-PLL','Location','northwest'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',white) %filling colour of legend
    set(h1,'Box','on') %Remove outer frame of legend
    ylim([-5 10])
    yticks(-5:1:10); % 设置Y轴刻度
    grid on
    box off
    title('KFs vs PLL')
    xlabel('Time step')
    ylabel('Phase[rad]')

end