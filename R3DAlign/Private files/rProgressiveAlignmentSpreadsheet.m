%The nucleotides specified by NTList1 of the structure in File1 have been
%aligned to the nucleotides specified by NTList2 of File2 and nucleotides
%NTList3 of File3 (separately).  AlignedNTList1 contains the indices of the
%nucleotides of File1 that are aligned with nucleotides in File2
%(AlignedNTList2).  AlignedNTList3 contains the indices of the nucleotides
%of File1 that are aligned with nucleotides in File 3 (AlignedNTList4).

function [] = rProgressiveAlignmentSpreadsheet(File1,NTList1,File2,NTList2,File3,NTList3,AlignedNTList1,AlignedNTList2,AlignedNTList3,AlignedNTList4,sheetname)

if ischar(File1),
  Filename = File1;
  File1 = zGetNTData(Filename,0);
end,
if ischar(File2),
  Filename = File2;
  File2 = zGetNTData(Filename,0);
end
if ischar(File3),
  Filename = File3;
  File3 = zGetNTData(Filename,0);
end

if nargin<11
  sheetname=[];
end

[i j k] = find(File1.Edge(NTList1,NTList1));      % find interactions in first file, k is interaction type

BPAI=zeros(length(i),2);                % will contain indices of nt's with interacting as specified
                                        % will be at most length(i) pairs
                                        % (actually much smaller)
ct=0;

%Determine which interactions are true basepairs.  To include all near
%basepairs, use alternate line
for m = 1:length(i)
%     if (abs(k(m)) <= 20 || abs(k(m)) >=24) && (abs(k(m)) <= 120 || abs(k(m)) >= 124)
      if (abs(k(m)) < 14)
        ct=ct+1;
        BPAI(ct,1)=NTList1(i(m));
        BPAI(ct,2)=NTList1(j(m)); 
      end
end
BPAI=BPAI(1:ct,:);         % Were only ct basepairs, not length(i)

%Do the same for second file
[i j k] = find(File2.Edge(NTList2,NTList2));
BPBI=zeros(length(i),2);
ct=0;
for m = 1:length(i)
%      if (abs(k(m)) <= 20 || abs(k(m)) >=24) && (abs(k(m)) <= 120 || abs(k(m)) >= 124)
    if (abs(k(m)) < 14)
        ct=ct+1;
        BPBI(ct,1)=NTList2(i(m));
        BPBI(ct,2)=NTList2(j(m));
    end
end
BPBI=BPBI(1:ct,:);

%Do the same for third file
[i j k] = find(File3.Edge(NTList3,NTList3));
BPCI=zeros(length(i),2);
ct=0;
for m = 1:length(i)
%      if (abs(k(m)) <= 20 || abs(k(m)) >=24) && (abs(k(m)) <= 120 || abs(k(m)) >= 124)
    if (abs(k(m)) < 14)
        ct=ct+1;
        BPCI(ct,1)=NTList3(i(m));
        BPCI(ct,2)=NTList3(j(m));
    end
end
BPCI=BPCI(1:ct,:);

BPList1=zeros(length(BPAI)+length(BPBI)+length(BPCI),2);
BPList2=zeros(length(BPAI)+length(BPBI)+length(BPCI),2);
BPList3=zeros(length(BPAI)+length(BPBI)+length(BPCI),2);
ct=0;

%BPList1 will contain the indices of the nucleotides making basepairs in A while
%BPList2 will contain the indices of the aligned nucleotides in B.  For
%example if nucleotides 100 and 300 A are making a basepair and are aligned to 101 and 302,
%respectively, then for some i BPList1(i,:)=[100 300] and BPList2(i,:)=[101 302].  If a nucleotide 
%is aligned to a gap, then a 0 is used to represent the gap.

for i=1:length(BPAI)        %Go through all basepairs in A
   BPList1(i,1)=BPAI(i,1);
   BPList1(i,2)=BPAI(i,2);
   ct=ct+1;
   p1=find(AlignedNTList1==BPAI(i,1));  %What 1st nt in bp is aligned to, if anything
   r1=find(AlignedNTList1==BPAI(i,2));  %What 2nd nt in bp is aligned to, if anything
   if isempty(p1)                % 1st nt isn't aligned to anything
      BPList2(i,1)=0;           
      if isempty(r1)             
         BPList2(i,2)=0;
      else
         BPList2(i,2)=AlignedNTList2(r1);
      end
   else                         %1st nt is aligned to a nt in B
      BPList2(i,1)=AlignedNTList2(p1); 
      if isempty(r1)
         BPList2(i,2)=0;
      else
         BPList2(i,2)=AlignedNTList2(r1); 
         q1=find(BPBI(:,1)==AlignedNTList2(p1));
         s1=find(BPBI(:,2)==AlignedNTList2(r1));
         t1=intersect(q1,s1);  %location of aligned bp in BPBI
         if ~isempty(t1)
            BPBI(t1,:)=[];   %remove it to avoid redundancy when looping through B basepairs
         end
      end
   end
   
   clear p1 r1 q1 s1 t1;
   
   p2=find(AlignedNTList3==BPAI(i,1));  %What 1st nt in bp is aligned to, if anything
   r2=find(AlignedNTList3==BPAI(i,2));  %What 2nd nt in bp is aligned to, if anything
   if isempty(p2)                % 1st nt isn't aligned to anything
      BPList3(i,1)=0;           
      if isempty(r2)             
         BPList3(i,2)=0;
      else
         BPList3(i,2)=AlignedNTList4(r2);
      end
   else                         %1st nt is aligned to a nt in B
      BPList3(i,1)=AlignedNTList4(p2); 
      if isempty(r2)
         BPList3(i,2)=0;
      else
         BPList3(i,2)=AlignedNTList4(r2); 
         q2=find(BPCI(:,1)==AlignedNTList4(p2));
         s2=find(BPCI(:,2)==AlignedNTList4(r2));
         t2=intersect(q2,s2);  %location of aligned bp in BPCI
         if ~isempty(t2)
            BPCI(t2,:)=[];   %remove it to avoid redundancy when looping through C basepairs
         end
      end
   end
end

ct=length(BPAI);

%Now go through basepairs of B that were not found to have a full match in
%A in the previous loop
for i=1:length(BPBI)
   ct=ct+1;
   BPList2(ct,1)=BPBI(i,1);
   BPList2(ct,2)=BPBI(i,2);
   p1=find(AlignedNTList2==BPBI(i,1));
   r1=find(AlignedNTList2==BPBI(i,2));
   if isempty(p1)
      BPList1(ct,1)=0;
      if isempty(r1)
         BPList1(ct,2)=0;
      else
         BPList1(ct,2)=AlignedNTList1(r1);
         y=find(AlignedNTList3==AlignedNTList1(r1));
         if isempty(y)
            BPList3(ct,2)=0;
         else
            BPList3(ct,2)=AlignedNTList4(y);
         end
      end
   else
      if isempty(r1)
         BPList1(ct,1)=AlignedNTList1(p1);
         BPList1(ct,2)=0;
         y=find(AlignedNTList3==AlignedNTList1(p1));
         if isempty(y)
            BPList3(ct,1)=0;
         else
            BPList3(ct,1)=AlignedNTList4(y);
         end
      else
         BPList1(ct,1)=AlignedNTList1(p1);
         BPList1(ct,2)=AlignedNTList1(r1); 
         y=find(AlignedNTList3==AlignedNTList1(p1));
         if isempty(y)
            BPList3(ct,1)=0;
         else
            BPList3(ct,1)=AlignedNTList4(y);
         end
         y=find(AlignedNTList3==AlignedNTList1(r1));
         if isempty(y)
            BPList3(ct,2)=0;        %changed from BPList3(ct,1)=0;
         else
            BPList3(ct,2)=AlignedNTList4(y);     %changed from BPList3(ct,1)=AlignedNTList4(y);
         end
         q1=find(BPAI(:,1)==AlignedNTList1(p1));
         s1=find(BPAI(:,2)==AlignedNTList1(r1));
         t1=intersect(q1,s1);
         if ~isempty(t1)
            BPAI(t1,:)=[];
         end
      end
   end
end

%Now go through basepairs of C that were not found to have a full match in
%A in the previous loop
for i=1:length(BPCI)
   ct=ct+1;
   BPList3(ct,1)=BPCI(i,1);
   BPList3(ct,2)=BPCI(i,2);
   p1=find(AlignedNTList4==BPCI(i,1));
   r1=find(AlignedNTList4==BPCI(i,2));
   if isempty(p1)
      BPList1(ct,1)=0;
      if isempty(r1)
         BPList1(ct,2)=0;
      else
         BPList1(ct,2)=AlignedNTList3(r1);
         y=find(AlignedNTList1==AlignedNTList3(r1));
         if isempty(y)
            BPList2(ct,2)=0;
         else
            BPList2(ct,2)=AlignedNTList2(y);
         end
      end
   else
      if isempty(r1)
         BPList1(ct,1)=AlignedNTList3(p1);  %changed this from AlignedNTList1(p1)
         BPList1(ct,2)=0;
         y=find(AlignedNTList1==AlignedNTList3(p1));
         if isempty(y)
            BPList2(ct,1)=0;
         else
            BPList2(ct,1)=AlignedNTList2(y);
         end
      else
         BPList1(ct,1)=AlignedNTList3(p1);      %....
         BPList1(ct,2)=AlignedNTList3(r1);      %....
         y=find(AlignedNTList1==AlignedNTList3(p1));
         if isempty(y)
            BPList2(ct,1)=0;
         else
            BPList2(ct,1)=AlignedNTList2(y);
         end
         y=find(AlignedNTList1==AlignedNTList3(r1));
         if isempty(y)
            BPList2(ct,2)=0;            %changed from BPList2(ct,1)=0;
         else
            BPList2(ct,2)=AlignedNTList2(y);  %changed from BPList2(ct,1)=AlignedNTList2(y);
         end
         q1=find(BPAI(:,1)==AlignedNTList3(p1));    %....
         s1=find(BPAI(:,2)==AlignedNTList3(r1));    %....
         t1=intersect(q1,s1);
         if ~isempty(t1)
            BPAI(t1,:)=[];
         end
      end
   end
end

%Add in aligned nucleotides not involved in any basepairs
for i = 1:length(AlignedNTList1)
   if ~any(BPList1(:,1)==AlignedNTList1(i)) || ~any(BPList2(:,1)==AlignedNTList2(i))
      ct=ct+1;
      BPList1(ct,1)=AlignedNTList1(i);
      BPList1(ct,2)=0;
      BPList2(ct,1)=AlignedNTList2(i);
      BPList2(ct,2)=0;
   end
   y=find(AlignedNTList3==AlignedNTList1(i));
   if isempty(y)
      BPList3(ct,1)=0;
      BPList3(ct,2)=0;
   else
      BPList3(ct,1)=AlignedNTList4(y);
      BPList3(ct,2)=0;
   end
end
for i = 1:length(AlignedNTList3)
   if ~any(BPList1(:,1)==AlignedNTList3(i)) || ~any(BPList3(:,1)==AlignedNTList4(i))
      ct=ct+1;
      BPList1(ct,1)=AlignedNTList3(i);
      BPList1(ct,2)=0;
      BPList3(ct,1)=AlignedNTList4(i);
      BPList3(ct,2)=0;
   end
   y=find(AlignedNTList1==AlignedNTList3(i));
   if isempty(y)
      BPList2(ct,1)=0;
      BPList2(ct,2)=0;
   else
      BPList2(ct,1)=AlignedNTList2(y);
      BPList2(ct,2)=0;
   end
end

%Add in unaligned nucleotides not involved in any basepairs
for i = 1:length(NTList1)
   if ~any(BPList1(:,1)==NTList1(i))
      ct=ct+1;
      BPList1(ct,1)=NTList1(i);
      BPList1(ct,2)=0;
      BPList2(ct,1)=0;
      BPList2(ct,2)=0;
      BPList3(ct,1)=0;
      BPList3(ct,2)=0;
   end
end
for i = 1:length(NTList2)
   if ~any(BPList2(:,1)==NTList2(i))
      ct=ct+1;
      BPList1(ct,1)=0;
      BPList1(ct,2)=0;
      BPList2(ct,1)=NTList2(i);
      BPList2(ct,2)=0;
      BPList3(ct,1)=0;
      BPList3(ct,2)=0;
   end
end
for i = 1:length(NTList3)
   if ~any(BPList3(:,1)==NTList3(i))
      ct=ct+1;
      BPList1(ct,1)=0;
      BPList1(ct,2)=0;
      BPList2(ct,1)=0;
      BPList2(ct,2)=0;
      BPList3(ct,1)=NTList3(i);
      BPList3(ct,2)=0;
   end
end

BPList1=BPList1(1:ct,:);
BPList2=BPList2(1:ct,:);
BPList3=BPList3(1:ct,:);

%Sort so nucleotides in first column are in order
[V,I]=sort(BPList1(:,1),'ascend');
BPList1=BPList1(I,:);
BPList2=BPList2(I,:);
BPList3=BPList3(I,:);

%At this point, all basepairs in B with the first nt aligned with a gap are
% found at the beginning of the list (since the first column of BPList1 is
% a 0).  These are moved to the locations so that the first column of
% BPList2 is ascending as well.
ct=1;
while BPList1(ct,1)==0 && BPList2(ct,1)>0
    ct = ct + 1;
end

num2move=ct-1;
N=length(BPList1(:,1));
for i=1:num2move
   tmpBP1=BPList1(1,:);
   tmpBP2=BPList2(1,:);
   tmpBP3=BPList3(1,:);
   P=find(BPList2(num2move+1:N,1)>tmpBP2(1,1),1,'first');
   if isempty(P)
      BPList1(1:N-1,:)=BPList1(2:N,:);
      BPList1(N,:)=tmpBP1;
      BPList2(1:N-1,:)=BPList2(2:N,:);
      BPList2(N,:)=tmpBP2;
      BPList3(1:N-1,:)=BPList3(2:N,:);
      BPList3(N,:)=tmpBP3;
   else
      BPList1(1:num2move+P-2,:)=BPList1(2:num2move+P-1,:);
      BPList1(num2move+P-1,:)=tmpBP1;
      BPList2(1:num2move+P-2,:)=BPList2(2:num2move+P-1,:);
      BPList2(num2move+P-1,:)=tmpBP2;
      BPList3(1:num2move+P-2,:)=BPList3(2:num2move+P-1,:);
      BPList3(num2move+P-1,:)=tmpBP3;
   end
   num2move=num2move-1;
end

%At this point, all basepairs in C with the first nt aligned with a gap are
% found at the beginning of the list (since the first column of BPList1 is
% a 0).  These are moved to the locations so that the first column of
% BPList3 is ascending as well.
ct=1;
while BPList1(ct,1)==0 && BPList3(ct,1)>0
    ct = ct + 1;
end

num2move=ct-1;
N=length(BPList1(:,1));
for i=1:num2move
   tmpBP1=BPList1(1,:);
   tmpBP2=BPList2(1,:);
   tmpBP3=BPList3(1,:);
   P=find(BPList3(num2move+1:N,1)>tmpBP3(1,1),1,'first');
   if isempty(P)
      BPList1(1:N-1,:)=BPList1(2:N,:);
      BPList1(N,:)=tmpBP1;
      BPList2(1:N-1,:)=BPList2(2:N,:);
      BPList2(N,:)=tmpBP2;
      BPList3(1:N-1,:)=BPList3(2:N,:);
      BPList3(N,:)=tmpBP3;
   else
      BPList1(1:num2move+P-2,:)=BPList1(2:num2move+P-1,:);
      BPList1(num2move+P-1,:)=tmpBP1;
      BPList2(1:num2move+P-2,:)=BPList2(2:num2move+P-1,:);
      BPList2(num2move+P-1,:)=tmpBP2;
      BPList3(1:num2move+P-2,:)=BPList3(2:num2move+P-1,:);
      BPList3(num2move+P-1,:)=tmpBP3;
   end
   num2move=num2move-1;
end


%At this point, all basepairs in B with the first nt aligned with a gap are
% found at the beginning of the list (since the first column of BPList1 is
% a 0).  These are moved to the locations so that the first column of
% BPList2 is ascending as well.
ct=1;
while BPList1(ct,1)==0 && BPList2(ct,1)>0
    ct = ct + 1;
end

num2move=ct-1;
N=length(BPList1(:,1));
for i=1:num2move
   tmpBP1=BPList1(1,:);
   tmpBP2=BPList2(1,:);
   tmpBP3=BPList3(1,:);
   P=find(BPList2(num2move+1:N,1)>tmpBP2(1,1),1,'first');
   if isempty(P)
      BPList1(1:N-1,:)=BPList1(2:N,:);
      BPList1(N,:)=tmpBP1;
      BPList2(1:N-1,:)=BPList2(2:N,:);
      BPList2(N,:)=tmpBP2;
      BPList3(1:N-1,:)=BPList3(2:N,:);
      BPList3(N,:)=tmpBP3;
   else
      BPList1(1:num2move+P-2,:)=BPList1(2:num2move+P-1,:);
      BPList1(num2move+P-1,:)=tmpBP1;
      BPList2(1:num2move+P-2,:)=BPList2(2:num2move+P-1,:);
      BPList2(num2move+P-1,:)=tmpBP2;
      BPList3(1:num2move+P-2,:)=BPList3(2:num2move+P-1,:);
      BPList3(num2move+P-1,:)=tmpBP3;
   end
   num2move=num2move-1;
end

%At this point, all basepairs in C with the first nt aligned with a gap are
% found at the beginning of the list (since the first column of BPList1 is
% a 0).  These are moved to the locations so that the first column of
% BPList3 is ascending as well.
ct=1;
while BPList1(ct,1)==0 && BPList3(ct,1)>0
    ct = ct + 1;
end

num2move=ct-1;
N=length(BPList1(:,1));
for i=1:num2move
   tmpBP1=BPList1(1,:);
   tmpBP2=BPList2(1,:);
   tmpBP3=BPList3(1,:);
   P=find(BPList3(num2move+1:N,1)>tmpBP3(1,1),1,'first');
   if isempty(P)
      BPList1(1:N-1,:)=BPList1(2:N,:);
      BPList1(N,:)=tmpBP1;
      BPList2(1:N-1,:)=BPList2(2:N,:);
      BPList2(N,:)=tmpBP2;
      BPList3(1:N-1,:)=BPList3(2:N,:);
      BPList3(N,:)=tmpBP3;
   else
      BPList1(1:num2move+P-2,:)=BPList1(2:num2move+P-1,:);
      BPList1(num2move+P-1,:)=tmpBP1;
      BPList2(1:num2move+P-2,:)=BPList2(2:num2move+P-1,:);
      BPList2(num2move+P-1,:)=tmpBP2;
      BPList3(1:num2move+P-2,:)=BPList3(2:num2move+P-1,:);
      BPList3(num2move+P-1,:)=tmpBP3;
   end
   num2move=num2move-1;
end

FinalListing=cell(length(BPList1(:,1)),8);

for i=1:length(FinalListing)
     if BPList1(i,1)==0
         a='---';
     else
         a = [File1.NT(BPList1(i,1)).Base File1.NT(BPList1(i,1)).Number];
     end
     if BPList1(i,1)==0 || BPList1(i,2)==0
        b=' ';
     else
         b = zEdgeText(File1.Edge(BPList1(i,1),BPList1(i,2)));
     end
     if BPList1(i,2)==0
         c='---';
     else
         c = [File1.NT(BPList1(i,2)).Base File1.NT(BPList1(i,2)).Number];
     end
     if BPList2(i,1)==0
         d='---';
     else
         d = [File2.NT(BPList2(i,1)).Base File2.NT(BPList2(i,1)).Number];
     end
     if BPList2(i,1)==0 || BPList2(i,2)==0
        e=' ';
     else
        e = zEdgeText(File2.Edge(BPList2(i,1),BPList2(i,2)));
     end
     if BPList2(i,2)==0
         f='---';
     else
         f = [File2.NT(BPList2(i,2)).Base File2.NT(BPList2(i,2)).Number];
     end
     if BPList3(i,1)==0
         g='---';
     else
         g = [File3.NT(BPList3(i,1)).Base File3.NT(BPList3(i,1)).Number];
     end
     if BPList3(i,1)==0 || BPList3(i,2)==0
        h=' ';
     else
        h = zEdgeText(File3.Edge(BPList3(i,1),BPList3(i,2)));
     end
     if BPList3(i,2)==0
         j='---';
     else
         j = [File3.NT(BPList3(i,2)).Base File3.NT(BPList3(i,2)).Number];
     end
 
     FinalListing{i,1}=a;
     FinalListing{i,2}=b;
     FinalListing{i,3}=c;
     FinalListing{i,4}=d;
     FinalListing{i,5}=e;
     FinalListing{i,6}=f;
     FinalListing{i,7}=g;
     FinalListing{i,8}=h;
     FinalListing{i,9}=j;
end


date = regexprep(datestr(now),':', '-');
ExcelName=['R3D Align ' File1.Filename ' ' File2.Filename ' ' File3.Filename ' ' date(1:17)];

FinalListing(2:end+1,:)=FinalListing(1:end,:);
FinalListing{1,1}=[];
FinalListing{1,2}=File1.Filename;
FinalListing{1,3}=[];
FinalListing{1,4}=[];
FinalListing{1,5}=File2.Filename;
FinalListing{1,6}=[];
FinalListing{1,7}=[];
FinalListing{1,8}=File3.Filename;
FinalListing{1,9}=[];
FinalListing{1,10}=[];

try

if isempty(sheetname)
   xlswrite([pwd filesep ExcelName],FinalListing)
else
   xlswrite([pwd filesep ExcelName],FinalListing,sheetname)
end


Excel = actxserver('Excel.Application');
Excel.Workbooks.Open([pwd '\' ExcelName]);

%Color basepairs according to geometric family
for j=[2 5 8]
   for i=1:length(FinalListing)
      if j==2
         Range = Excel.Range(['B' num2str(i)]);
      elseif j==5
         Range = Excel.Range(['E' num2str(i)]);
      else
         Range = Excel.Range(['H' num2str(i)]);
      end
      switch lower(strtrim(FinalListing{i,j}))
         case 'cww'
            Range.Interior.ColorIndex = 4;
         case {'tsh','ths'}
            Range.Interior.ColorIndex = 38;
         case {'cws','csw'}
            Range.Interior.ColorIndex = 33;
         case {'tws','tsw'}
            Range.Interior.ColorIndex = 41;     
         case 'tss'
            Range.Interior.ColorIndex = 22;
         case {'thw','twh'}
            Range.Interior.ColorIndex = 45; 
         case 'css'
            Range.Interior.ColorIndex = 36;
         case {'chs','csh'}
            Range.Interior.ColorIndex = 7; 
         case {'chw','cwh'}
            Range.Interior.ColorIndex = 40;
         case 'thh'
            Range.Interior.ColorIndex = 12;
         case 'tww'
            Range.Interior.ColorIndex = 6;   
         case 'chh'
            Range.Interior.ColorIndex = 9;
      end
   end
end

%Color near-basepairs according to geometric family
for j=[2 5 8]
   for i=1:length(FinalListing)
      if j==2
         Range = Excel.Range(['B' num2str(i)]);
      elseif j==5
         Range = Excel.Range(['E' num2str(i)]);
      else
         Range = Excel.Range(['H' num2str(i)]);
      end
      switch lower(strtrim(FinalListing{i,j}))
         case 'ncww'
            Range.Interior.ColorIndex = 4;
         case {'ntsh','nths'}
            Range.Interior.ColorIndex = 38;
         case {'ncws','ncsw'}
            Range.Interior.ColorIndex = 33;
         case {'ntws','ntsw'}
            Range.Interior.ColorIndex = 41;     
         case 'ntss'
            Range.Interior.ColorIndex = 22;
         case {'nthw','ntwh'}
            Range.Interior.ColorIndex = 45; 
         case 'ncss'
            Range.Interior.ColorIndex = 36;
         case {'nchs','ncsh'}
            Range.Interior.ColorIndex = 7; 
         case {'nchw','ncwh'}
            Range.Interior.ColorIndex = 40;
         case 'nthh'
            Range.Interior.ColorIndex = 12;
         case 'ntww'
            Range.Interior.ColorIndex = 6;   
         case 'nchh'
            Range.Interior.ColorIndex = 9;
      end
   end
end

catch eObj1
   eObj1.identifier
   Excel.Workbooks.Close;
   Excel.Quit;
   Excel.delete;
   throw(eObj1)
end

Excel.Quit;
Excel.delete;