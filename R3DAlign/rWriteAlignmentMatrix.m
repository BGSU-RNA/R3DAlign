function [MA] = rWriteAlignmentMatrix(File1,NTIndices1,File2,NTIndices2,AlignedNTIndices1,AlignedNTIndices2,NTList,filename)

if ischar(File1),
  Filename = File1;
  File1 = zAddNTData(Filename);
end
if ischar(File2),
  Filename = File2;
  File2 = zAddNTData(Filename);
end

AlignmentName=[filename '.txt'];
fidOUT = fopen(AlignmentName,'w+');

numAl=length(AlignedNTIndices1);
Alignment=zeros(length(NTIndices1)+length(NTIndices2)-numAl,2);

for i=1:length(NTIndices1)
   t=find(AlignedNTIndices1==NTIndices1(i));
   if isempty(t)
      Alignment(i,:)=[NTIndices1(i) 0];
   else
      Alignment(i,:)=[NTIndices1(i) AlignedNTIndices2(t)];
   end
end

ct=length(NTIndices1);

for i=NTIndices2
   if ~any(AlignedNTIndices2==i);
       s=find(Alignment(:,2)>i,1,'first');
       if isempty(s)
          ct=ct+1;
          Alignment(ct,:)=[0 i];
       else          
          Alignment(s+1:ct+1,:)=Alignment(s:ct,:);
          Alignment(s,:)=[0 i];
          ct=ct+1;
       end
   end
end

L = length(Alignment);

line=cell(1,2);
%  if ischar(NTList{1}) && ~isequal(lower(NTList{1}),'all')
%     line{1} =[File1.Filename ' ' NTList{1} ' '];
%   %  fprintf(fidOUT,'> %s %s\r\n',File1.Filename,NTList{1});
%  elseif iscell(NTList{1})
%     line{1} =[File1.Filename ' ' NTList{1}{1} ' '];
%  %   fprintf(fidOUT,'> %s %s\r\n',File1.Filename,NTList{1}{1});
%  else
%     line{1} =[File1.Filename ' '];
%   %  fprintf(fidOUT,'> %s\r\n',File1.Filename);
%  end

%  if ischar(NTList{2}) && ~isequal(lower(NTList{2}),'all')
%     line{2} =[File2.Filename ' ' NTList{2} ' '];
%  %   fprintf(fidOUT,'> %s %s\r\n',File2.Filename,NTList{2});
%  elseif iscell(NTList{2})
%     line{2} =[File2.Filename ' ' NTList{2}{1} ' '];
%   %  fprintf(fidOUT,'> %s %s\r\n',File2.Filename,NTList{2}{1});
%  else
%     line{2} =[File2.Filename ' '];
%    % fprintf(fidOUT,'> %s\r\n',File2.Filename);
%  end

% The above was commented out and replaced with the two lines below because the 
% nt numbers were being written as 200,201,202... instead of 200:300 (for example).
% To correct this, the text entered would have to be saved before being parsed.

  line{1} =[File1.Filename ' '];
  line{2} =[File2.Filename ' '];

  fprintf(fidOUT,'\r\n');
  a=length(line{1});
  b=length(line{2});

  [m1 n1] = min([a b]);

  for i = 1:abs(a-b)
     line{n1} = [line{n1} ' '];
  end
  
  c1 = 1;
  C1 = 1;
  c2 = 1;
  C2 = 1;
  while c1 <= (80-max(a,b)) && C1 <=L
    if Alignment(C1,1)==0
       line{1} = [line{1} '-'];
    else
       line{1} = [line{1} File1.NT(Alignment(C1,1)).Base];
    end
    if Alignment(C1,2)==0
       line{2} = [line{2} '-'];
    else
       line{2} = [line{2} File2.NT(Alignment(C1,2)).Base];
    end
    c1 = c1 + 1;
    C1 = C1 + 1;
    c2 = c2 + 1;
    C2 = C2 + 1;
  end
  fprintf(fidOUT,'%s\r\n',line{1});
  fprintf(fidOUT,'%s\r\n',line{2});
  fprintf(fidOUT,'\r\n');

  while C1<=L
    c1=1;
    while c1<=80 && C1<=L
       if Alignment(C1,1)==0
          fprintf(fidOUT,'-');
       else
          fprintf(fidOUT,'%s', File1.NT(Alignment(C1,1)).Base);
       end
       c1 = c1 + 1;
       C1 = C1 + 1;
    end
    fprintf(fidOUT,'\r\n');

    c2=1;
    while c2<=80 && C2<=L
       if Alignment(C2,2)==0
          fprintf(fidOUT,'-');
       else
          fprintf(fidOUT,'%s', File2.NT(Alignment(C2,2)).Base);
       end
       c2 = c2 + 1;
       C2 = C2 + 1;
    end
    fprintf(fidOUT,'\r\n');
    fprintf(fidOUT,'\r\n');
  end
    
fclose(fidOUT);