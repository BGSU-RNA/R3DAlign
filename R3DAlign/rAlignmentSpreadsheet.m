function [] = rAlignmentSpreadsheet(File1,NTList1,File2,NTList2,AlignedNTList1,AlignedNTList2,sheetname)

if ischar(File1),
  Filename = File1;
  File1 = zGetNTData(Filename,0);
end
if ischar(File2),
  Filename = File2;
  File2 = zGetNTData(Filename,0);
end

if nargin<7
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

BPList1=zeros(length(BPAI)+length(BPBI),2);
BPList2=zeros(length(BPAI)+length(BPBI),2);
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
   p=find(AlignedNTList1==BPAI(i,1));  %What 1st nt in bp is aligned to, if anything
   r=find(AlignedNTList1==BPAI(i,2));  %What 2nd nt in bp is aligned to, if anything
   if isempty(p)                % 1st nt isn't aligned to anything
      BPList2(i,1)=0;           
      if isempty(r)             
         BPList2(i,2)=0;
      else
         BPList2(i,2)=AlignedNTList2(r);
      end
   else                         %1st nt is aligned to a nt in B
      BPList2(i,1)=AlignedNTList2(p); 
      if isempty(r)
         BPList2(i,2)=0;
      else
         BPList2(i,2)=AlignedNTList2(r); 
         q=find(BPBI(:,1)==AlignedNTList2(p));
         s=find(BPBI(:,2)==AlignedNTList2(r));
         t=intersect(q,s);  %location of aligned bp in BPBI
         if ~isempty(t)
            BPBI(t,:)=[];   %remove it to avoid redundancy when looping through B basepairs
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
   p=find(AlignedNTList2==BPBI(i,1));
   r=find(AlignedNTList2==BPBI(i,2));
   if isempty(p)
      BPList1(ct,1)=0;
      if isempty(r)
         BPList1(ct,2)=0;
      else
         BPList1(ct,2)=AlignedNTList1(r);
      end
   else
      if isempty(r)
         BPList1(ct,1)=AlignedNTList1(p);
         BPList1(ct,2)=0;
      else
         BPList1(ct,1)=AlignedNTList1(p);
         BPList1(ct,2)=AlignedNTList1(r); 
         q=find(BPAI(:,1)==AlignedNTList1(p));
         s=find(BPAI(:,2)==AlignedNTList1(r));
         t=intersect(q,s);
         if ~isempty(t)
            BPAI(t,:)=[];
         end
      end
   end
end

%Add it aligned nucleotides not involved in any basepairs
for i = 1:length(AlignedNTList1)
   if ~any(BPList1(:,1)==AlignedNTList1(i)) || ~any(BPList2(:,1)==AlignedNTList2(i))
      ct=ct+1;
      BPList1(ct,1)=AlignedNTList1(i);
      BPList1(ct,2)=0;
      BPList2(ct,1)=AlignedNTList2(i);
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
   end
end
for i = 1:length(NTList2)
   if ~any(BPList2(:,1)==NTList2(i))
      ct=ct+1;
      BPList1(ct,1)=0;
      BPList1(ct,2)=0;
      BPList2(ct,1)=NTList2(i);
      BPList2(ct,2)=0;
   end
end
BPList1=BPList1(1:ct,:);
BPList2=BPList2(1:ct,:);

%Sort so nucleotides in first column are in order
[V,I]=sort(BPList1(:,1),'ascend');
BPList1=BPList1(I,:);
BPList2=BPList2(I,:);


%At this point, all basepairs in B with the first nt aligned with a gap are
% found at the beginning of the list (since the first column of BPList1 is
% a 0).  These are moved to the locations so that the first column of
% BPList2 is ascending as well.
num2move=length(find(BPList1(:,1)==0));
N=length(BPList1(:,1));
for i=1:num2move
   tmpBP1=BPList1(1,:);
   tmpBP2=BPList2(1,:);
   P=find(BPList2(num2move+1:N,1)>tmpBP2(1,1),1,'first');
   if isempty(P)
      BPList1(1:N-1,:)=BPList1(2:N,:);
      BPList1(N,:)=tmpBP1;
      BPList2(1:N-1,:)=BPList2(2:N,:);
      BPList2(N,:)=tmpBP2;
   else
      BPList1(1:num2move+P-2,:)=BPList1(2:num2move+P-1,:);
      BPList1(num2move+P-1,:)=tmpBP1;
      BPList2(1:num2move+P-2,:)=BPList2(2:num2move+P-1,:);
      BPList2(num2move+P-1,:)=tmpBP2;
   end
   num2move=num2move-1;
end

FinalListing=cell(length(BPList1(:,1)),6);

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
         
     FinalListing{i,1}=a;
     FinalListing{i,2}=b;
     FinalListing{i,3}=c;
     FinalListing{i,4}=d;
     FinalListing{i,5}=e;
     FinalListing{i,6}=f;
end

date = regexprep(datestr(now),':', '-');
ExcelName=['R3D Align ' File1.Filename ' ' File2.Filename ' ' date(1:17)];

FinalListing(2:end+1,:)=FinalListing(1:end,:);
FinalListing{1,1}=[];
FinalListing{1,2}=File1.Filename;
FinalListing{1,3}=[];
FinalListing{1,4}=[];
FinalListing{1,5}=File2.Filename;
FinalListing{1,6}=[];

try
if isempty(sheetname)
   xlswrite([pwd filesep ExcelName],FinalListing)
else
   xlswrite([pwd filesep ExcelName],FinalListing,sheetname)
end

Excel = actxserver('Excel.Application');
Excel.Workbooks.Open([pwd '\' ExcelName]);

%Color basepairs according to geometric family
for j=[2 5]
   for i=1:length(FinalListing)
      if j==2
         Range = Excel.Range(['B' num2str(i)]);
      else
         Range = Excel.Range(['E' num2str(i)]);
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

catch eObj1
   eObj1.identifier
   Excel.Workbooks.Close;
   Excel.Quit;
   Excel.delete;
   throw(eObj1)
end

Excel.Quit;
Excel.delete;