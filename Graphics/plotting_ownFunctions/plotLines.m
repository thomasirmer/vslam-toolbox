function [] = plotLines(image, lines, hFigureImage)
    
    figure(hFigureImage);
    
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