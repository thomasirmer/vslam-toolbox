function [ h ] = plotHmgLine(e, style)
    % homogenous line coordinates
    l = e(1);
    m = e(2);
    n = e(3);
    
    % homogenous line
    x = 0:600;
    y = - (n + l * x) / m;
    
    % plotting
    h = plot(x, y, style);
end