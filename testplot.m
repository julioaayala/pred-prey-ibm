function testplot
    x = 1:200;
    Y = zeros(2,200);
    Y2 = zeros(2,200);
    Y(1,:) = log(x);
    Y(2,:) = log(x+10);
    Y2(1,:) = log(x.^-1);
    Y2(2,:) = log((x+10).^-1);
    figure;
    legends = union ("Pop" + string(1:2), ...
                     "Res" + string(1:2));
    h(1) = animatedline('Color','#0072BD');
    h(2) = animatedline('Color','#D95319');
    h(3) = animatedline('Color','#EDB120');
    h(4) = animatedline('Color','#7E2F8E');
    xlabel('t')
    ylabel('abun')
    legend(legends);
    for i=1:length(x)
        addpoints(h(1),x(i),Y(1,i));
        addpoints(h(2),x(i),Y2(1,i));
        addpoints(h(3),x(i),Y(2,i));
        addpoints(h(4),x(i),Y2(2,i));
        drawnow;
    end

end