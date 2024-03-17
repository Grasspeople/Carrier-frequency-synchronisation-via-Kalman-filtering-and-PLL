function draw_UKF(Nsteps,x_u_series,x_truth)
    orange = [1 0.34 0.20]; 
    lightgrey = [0.94 0.94 0.94]; % color definition
    blue = [0.21 0.35 1]; 
    white=[1 1 1];
    xlabel('steps')
    ylabel('phase [rad]')
    %%    
    figure(1)

    plot((1:Nsteps),x_truth(1,1:Nsteps),'.-','Color',blue,'LineWidth',3)
    hold on     
    plot((1:Nsteps),x_u_series(1,1:Nsteps),'.-','Color',orange,'LineWidth',2)
    h1 = legend('Truth','UKF','Location','southeast'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',white) %filling colour of legend
    set(h1,'Box','on') %Remove outer frame of legend
    grid on
    box off
    title('UKF')
    xlabel('Time step')
    ylabel('Phase[rad]')




end