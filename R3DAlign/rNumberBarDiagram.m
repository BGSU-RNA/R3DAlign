% View(8) = 1 means to space the chains and gaps in the normal way
% View(9) = 1 means to display the numbers; 0 means just to calculate

function [A,mA] = rNumberBarDiagram(File,NTList,View,Thickness,r,numLoc)

N = length(NTList);

% ------------------------------------- Determine where to plot each NT

A = zeros(1,N);      
A(1) = 1;

if View(8) == 1,
  spaces = [1 4 8 12];
else
  spaces = [1 0 0 12];
end

nc = 1;                               
for t = 1:(N-1),
  A(t+1) = A(t) + spaces(1);
  if File.Covalent(NTList(t),NTList(t+1)) == 0,
    A(t+1) = A(t+1) + spaces(2);
  end
  if (File.NT(NTList(t)).Chain ~= File.NT(NTList(t+1)).Chain),
    A(t+1) = A(t+1) + spaces(3);
    nc = [nc t+1];
  end
end
A(end+1) = A(end) + spaces(4);
mA = max(A);

barlength=20;
A = A * barlength / mA;

if View(9) == 1,
   plot([0 barlength],[r r], 'k');
   hold on;
   
   for n = 1:length(nc),
      theta = A(nc(n));
      if isequal(numLoc,'above')
         text(theta, r+.02, ['Chain ' File.NT(NTList(nc(n))).Chain], 'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 6);
      else
         text(theta, r-.02, ['Chain ' File.NT(NTList(nc(n))).Chain], 'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 6);
      end
   end
end

   if N > 1000,
      step  = 50;
      sstep = 10;
   elseif N > 500,
      step = 20;
      sstep = 5;
   elseif N > 300,
      step = 10;
      sstep = 2;
   elseif N > 100,
      step = 5;
      sstep = 1;
   else
      step = 2;
      sstep = 2;
   end

   kk = 1;

   flag = 0;
   for k = 1:N,
      kkk = str2num(File.NT(NTList(k)).Number);
      if ~isempty(kkk),                               % it's really a number
         kk = kkk;
         flag = 0;                                     % OK to use next w/ letter
      else                                            % it has a letter in it
         kkk = str2num(File.NT(NTList(k)).Number(1:(end-1)));  % omit last character
         if ~isempty(kkk),
            kk = kkk;                                   % use this for display
            flag = 1;                                   % but only once
         end
      end
      theta = A(k);

      if isequal(numLoc,'above')
         if (mod(kk,sstep) == 0) && (mod(kk,step) ~= 0) && (Thickness < 1) && (flag == 0),
            plot([theta theta], [r, r+.01],'k');
         end
         if mod(kk,step) == 0 && (flag == 0),
            plot([theta theta], [r r+.02],'k');
            text(theta, r+.02, File.NT(NTList(k)).Number,'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
         end
      elseif isequal(numLoc,'below')
         if (mod(kk,sstep) == 0) && (mod(kk,step) ~= 0) && (Thickness < 1) && (flag == 0),
            plot([theta theta], [r-.01, r],'k');
         end
         if mod(kk,step) == 0 && (flag == 0),
            plot([theta theta], [r-.02 r],'k');
            text(theta, r-.02, File.NT(NTList(k)).Number, 'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
         end
      end
   end
end