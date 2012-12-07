function [] = rWriteAlignmentFasta(File1,NTIndices1,File2,NTIndices2,AlignedNTIndices1,AlignedNTIndices2,NTList,filename)

if nargin < 8
    filename = 'Alignment'; 
end

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
for k=1:2,
  if k==1
    File=File1;
  else
    File=File2;
  end
%  if ischar(NTList{k}) && ~isequal(lower(NTList{k}),'all')
%     fprintf(fidOUT,'> %s %s\r\n',File.Filename,NTList{k});
%  elseif iscell(NTList{k})
%     fprintf(fidOUT,'> %s %s\r\n',File.Filename,NTList{k}{1});
%  else
%     fprintf(fidOUT,'> %s\r\n',File.Filename);
%  end
 
% The above was commented out and replaced with the line below because the 
% nt numbers were being written as 200,201,202... instead of 200:300 (for example).
% To correct this, the text entered would have to be saved before being parsed.

 fprintf(fidOUT,'> %s\r\n',File.Filename);

 C = 1;
  while C <= L
    c=1;
    while c<=80 && C<=L
       if Alignment(C,k)==0
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