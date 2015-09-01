function [] = plotSegLine(coord, style)
    % line segment coordinates
    line_seg_x = [coord(1), coord(3)];
    line_seg_y = [coord(2), coord(4)];
    
    plot(line_seg_x, line_seg_y, style);
end