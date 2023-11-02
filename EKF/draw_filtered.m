function draw_filtered(Nsteps,y_measure,x_truth,x_u_series)
    orange = [1 0.34 0.20]; 
    lightgrey = [0.94 0.94 0.94]; % color definition
    blue = [0.21 0.35 1]; 
    xlabel('steps')
    ylabel('phase [rad]')
    %%    
    figure(1)
%     subplot(1, 2, 1);
%    plot((1:Nsteps),y_measure(1,:),'.-','Color',orange)
 %   hold on
    plot((1:Nsteps),x_truth(1,:),'.-','Color',blue)
    hold on     
plot((1:Nsteps),x_u_series(1,:),'.-g')
     title('Truth vs EKF')
    legend('truth','EKF','Location','northwest'); 
    xlabel('Nsteps')
    ylabel('phase[rad]')
    %%
%     subplot(1, 2, 2); 
%     plot((1:Nsteps),y_measure(2,:),'.-','Color',orange)
%     hold on
    %plot((1:Nsteps),x_truth(2,:),'.-','Color',blue)
    %hold on
%     title('Imaginary part')
%     h1 = legend('Measurements','Truth','KF','Location','northwest'); 
%     legend('Measurement','Location','northwest');
%     xlabel('Nsteps')
%     ylabel('phase[rad]')
end