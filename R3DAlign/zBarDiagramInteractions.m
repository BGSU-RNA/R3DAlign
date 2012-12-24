% zBarDiagramInteractions(File,NTList,A,View,Location) draws circular arcs representing pairwise interactions, with nucleotides spaced along a line as given by A

function [m] = zBarDiagramInteractions(File,NTList,A,View,Location)

if length(NTList) > 1000,
  Thickness = 0.2;
else
  Thickness = 1;
end

E  = fix(abs(File.Edge(NTList,NTList)));
C  = File.Crossing(NTList,NTList);

Color = zKeepZeros(zSparseValues(E,1,1),C);      % nested cWW
Color = Color + 2*zKeepZeros(zSparseRange(E,2,13,1),C); % nested non-cWW
Color = Color + 3*(E==1).*(C>0);                 % non-nested cWW
Color = Color + 4*zSparseRange(E,2,13,1).*(C>0); % non-nested non-cWW
Color = Color + 5*zSparseRange(E,21,24,1);       % stacks

clear E C

[i,j,c] = find(triu(Color));            % keep one of each interaction

BP = abs(File.BasePhosphate(NTList,NTList));           % 
BP = zSparseRange(BP,0,99,1);           % extract BPh interactions
[ii,jj,cc] = find(6*BP);                % base-phosphate interactions

k = find(ii ~= jj);                     % eliminate self BPh interactions

i = [i; ii(k)];                         % append base-phosphate interactions
j = [j; jj(k)];
c = [c; cc(k)];

BR = abs(File.BaseRibose(NTList,NTList));              % 
BR = zSparseRange(BR,0,99,1);           % extract BR interactions
[ii,jj,cc] = find(7*BR);                % base-ribose interactions

k = find(ii ~= jj);                     % eliminate self BR interactions

i = [i; ii(k)];                         % append base-ribose interactions
j = [j; jj(k)];
c = [c; cc(k)];

% ---------------------------------------- Draw the interactions in color

[y,k] = sort(-abs(c));               % sort by decreasing color, for overlap

i = i(k);
j = j(k);
c = c(k);

color = 'bcrgymo';

switch Location
  case 'above' 
    x = A;
    y = 0.02 * ones(size(x));
    a = 0:0.01:pi;
  case 'below'
    x = A;
    y = -0.22 * ones(size(x));
    a = pi:0.01:2*pi;
end

m = [0 0 0 0];

% interaction i(k),j(k) should be plotted if i(k) and j(k) are in NTList
% A, x, and y map NTList to coordinates

for k = 1:length(i),
  if View(c(k)) > 0,
    mm = zBarDiagramArc([x(i(k)) x(j(k))], [y(i(k)) y(j(k))], a, c(k), Thickness);
    m(1) = min(m(1),mm(1));
    m(2) = max(m(2),mm(2));
    m(3) = min(m(3),mm(3));
    m(4) = max(m(4),mm(4));

    hold on
  end
end
