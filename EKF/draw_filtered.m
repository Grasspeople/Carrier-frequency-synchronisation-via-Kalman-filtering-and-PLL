function draw_filtered(Nsteps,x_u_series,x_truth)
    orange = [1 0.34 0.20]; 
    lightgrey = [0.94 0.94 0.94]; % color definition
    blue = [0.21 0.35 1]; 
    xlabel('steps')
    ylabel('phase [rad]')
    %%    
    figure(1)
    plot((1:Nsteps),x_truth(1,:),'.-','Color',blue)
    hold on     
    plot((1:Nsteps),x_u_series(1,:),'.-','Color',orange)
    title('Truth vs EKF')
    legend('truth','EKF','Location','northwest'); 
    xlabel('Nsteps')
    ylabel('phase[rad]')

end