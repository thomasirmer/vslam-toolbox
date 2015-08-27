function [] = plotHmgLine(meas, exp)
    % homogenous line coordinates
    l = exp.e(1);
    m = exp.e(2);
    n = exp.e(3);
    
    % line segment coordinates
    line_seg_x = [meas.y(1), meas.y(3)];
    line_seg_y = [meas.y(2), meas.y(4)];
    
    % homogenous line
    x = 0:600;
    y = - (n + l * x) / m;
    
    hold on
    plot(line_seg_x, line_seg_y, 'blue--*');
    plot(x, y, 'green--');
end