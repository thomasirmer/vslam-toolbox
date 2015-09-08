function [] = plotLines(image, lines, hFigureImage)
    
    figure(hFigureImage);
    
    screens = get(0,'MonitorPositions');
    left    = screens(2,1) + 50;
    bottom  = 100;
    width   = screens(2,3) - 100;
    height  = screens(2,4) - 200;
    
    imshow(image);
    
    axis on
    
    hold on;
    
    x = [lines(:,1)' ; lines(:,3)'];
    y = [lines(:,2)' ; lines(:,4)'];
    
    for i = 1:length(x)
        line([x(1,i), x(2,i)], [y(1,i), y(2,i)], 'Color', 'red');
        text((x(1,i) + x(2,i)) / 2, (y(1,i) + y(2,i)) / 2, num2str(i), 'Color', 'red');
    end
    
    hold off;
end