% zBarDiagramArc([x(i(k)) x(j(k))], [y(i(k)) y(j(k))], c(k), Thickness);

function [m] = zBarDiagramArc(x, y, a, c, Thickness);

center = [(x(1)+x(2))/2 (y(1)+y(2))/2];

radius = norm([x(1) y(1)] - [x(2) y(2)])/2;

switch c
case 1,
  h = [0 0 1];                          % blue
case 2,
  h = [0 1 1];                          % cyan
case 3,
  h = [1 0 0];                          % red
case 4,
  h = [0 1 0];                          % green
case 5,
  h = [245 245 15]/255;                 % dark yellow
case 6,
  h = [0.8 0 1];                        % magenta
case 7,
  h = [255 165 0]/255;                  % orange
otherwise,
  h = [0 0 0];                          % black
end

s = center(1) + radius*cos(a);
t = center(2) + 0.05*radius*sin(a);

m(1) = min(s);
m(2) = max(s);
m(3) = min(t);
m(4) = max(t);

plot(s, t, 'Color', h, 'LineWidth', Thickness);
