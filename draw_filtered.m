function draw_filtered(Nsteps,y_measure,x_truth,x_u_series)
    orange = [1 0.34 0.20]; 
    lightgrey = [0.94 0.94 0.94]; % color definition
    blue = [0.21 0.35 1]; 
    xlabel('steps')
    ylabel('phase [rad]')
    
    figure(1)
    plot((1:Nsteps),y_measure(1,:),'.-','Color',orange)
    hold on
    plot((1:Nsteps),x_truth(1,:),'.-','Color',blue)
    hold on
    plot((1:Nsteps),x_u_series(1,:),'.-g')

    h1 = legend('Measurements','Truth','KF','Location','northwest'); 
    set(gca,'linewidth',0.4); % thickness of grid
    set(gca,'GridLineStyle','-.');% type of grid
    set(gca,'GridAlpha',0.4); % dark of grid
    set(h1,'Color',lightgrey) %filling colour of legend
    set(h1,'Box','off') %Remove outer frame of legend
    grid on
    box off
    title('KF (green), Measurements (red), True (blue)')
end