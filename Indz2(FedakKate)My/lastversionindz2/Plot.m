function [ ] = Plot( params )
    figure('Position',[100 100 1000 300])
    subplot(1,2,1)
    plot(params.t,params.y(1,:))
    grid on
    title('y1, initial parameters')
    ylabel('y1')
    xlabel('time, t')
    subplot(1,2,2)
    plot(params.t,params.y(2,:))
    grid on
    title('y2, initial parameters')
    ylabel('y2')
    xlabel('time, t')
end

