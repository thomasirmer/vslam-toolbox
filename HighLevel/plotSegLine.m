function [] = plotHmgLine(coord)
    [hmgCoords, ~] = seg2hmgLin(coord(:,1));
    hmgCoords = hmgCoords ./ hmgCoords(1);
    
    % homogenous line coordinates
    l = hmgCoords(1);
    m = hmgCoords(2);
    n = hmgCoords(3);
    
    % line segment coordinates
    line_seg_x = [coord(1), coord(3)];
    line_seg_y = [coord(2), coord(4)];
    
    % homogenous line
    x = 0:600;
    y = - (n + l * x) / m;
    
    hold on
    plot(line_seg_x, line_seg_y, 'magenta--*');
    plot(x, y, 'magenta--');
end