function [] = rWriteAlignmentFasta(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,filename)

if ischar(File1),
  Filename = File1;
  File1 = zAddNTData(Filename);
end
if ischar(File2),
  Filename = File2;
  File2 = zAddNTData(Filename);
end

FastaName=[filename '.fasta'];
fidOUT = fopen(FastaName,'w+');

numAl=length(AlignedIndices1);
Alignment=zeros(length(Indices1)+length(Indices2)-numAl,2)-99999;

for i=1:length(Indices1)
   t=find(AlignedIndices1==Indices1(i));
   if isempty(t)
      Alignment(i,:)=[Indices1(i) -99999];
   else
      Alignment(i,:)=[Indices1(i) AlignedIndices2(t)];
   end
end

ct=length(Indices1);

for i=Indices2
   if ~any(AlignedIndices2==i);
       Diff=abs(i-Alignment(:,2));
       if ~isempty(find(Diff==0))
           disp('Problem')
           return;
       end
       [V loc]=min(Diff);
       if i-Alignment(loc,2) > 0
           Alignment(loc+2:ct+2,:)=Alignment(loc+1:ct+1,:);
           Alignment(loc+1,:)=[-99999 i];
           ct=ct+1;
       else
           Alignment(loc+1:ct+1,:)=Alignment(loc:ct,:);
           Alignment(loc,:)=[-99999 i];
           ct=ct+1;
       end
%        s=find(Alignment(:,2)>i,1,'first');
%        if isempty(s)
%           ct=ct+1;
%           Alignment(ct,:)=[0 i];
%        else          
%           Alignment(s+1:ct+1,:)=Alignment(s:ct,:);
%           Alignment(s,:)=[0 i];
%           ct=ct+1;
%        end
   end
end

L = length(Alignment);

for k=1:2,
  if k==1
    File=File1;
  else
    File=File2;
  end
%   if ischar(NTList{k}) && ~isequal(lower(NTList{k}),'all')
%      fprintf(fidOUT,'> %s %s\r\n',File(k).Filename,NTList{k});
%   elseif iscell(NTList{k})
%      fprintf(fidOUT,'> %s %s\r\n',File(k).Filename,NTList{k}{1});
%   else
%      fprintf(fidOUT,'> %s\r\n',File(k).Filename);
%   end

% The above was commented out and replaced with the line below because the 
% nt numbers were being written as 200,201,202... instead of 200:300 (for example).
% To correct this, the text entered would have to be saved before being parsed.

 fprintf(fidOUT,'> %s\r\n',File.Filename);
 
 C = 1;
  while C <= L
    c=1;
    while c<=80 && C<=L
       if Alignment(C,k)==-99999
          fprintf(fidOUT,'-');
       else
          fprintf(fidOUT,'%s', File.NT(Alignment(C,k)).Base);
       end
       c = c + 1;
       C = C + 1;
    end
    fprintf(fidOUT,'\r\n');
  end
end

fclose(fidOUT);