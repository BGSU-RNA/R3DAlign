% rAugmentAlignment(align1,align2) fills in a near match for non-aligned nucleotides
% This is in the seed alignment.
% For each nt in A aligned with a gap, the nearest aligned nt in B is found
% to use as the center of the band for that nt

% New in V2
% 1) Returns D so that R3DAlign can use D to set the anchors
% 2) More sophisticated approach to Augmenting alignment (instead of just
% finding nearest aligned, it "fills in" the correspondences on a one to one
% basis (as much as possible)).  For example if the seed alignment is:
% 5---6789
% 5678---9
% the resulting augmented alignment will be
% 56789
% 56789
% instead of
% 56789
% 59999
% like is was previously.

function [align1, align2, D] = rAugmentAlignment(Indices1,Indices2,align1,align2,D)

for i = 1:length(Indices1),
  if ~any(align1==i)
    align1 = [align1 i]; %#ok<AGROW>
    align2 = [align2 0]; %#ok<AGROW>
    D = [D 10]; %#ok<AGROW> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
end
[S1 IX] = sort(align1);
align1 = S1; %#ok<NASGU>
align2 = align2(IX);
D = D(IX);%%%%%%%%%%%%%%%%%%%

GappedAs=find(align2==0);

if sum(align2) == 0  %This happens when seed alignment did not align any
   align2=align1;
   return;
end

Diffs=[1 diff(GappedAs)];
currRow=1;
currCol=0;
for i=1:length(Diffs)
   if Diffs(i)~=1
      currRow=currRow+1;
      currCol=1;
   else
      currCol=currCol+1;
   end
   Mat(currRow,currCol)=GappedAs(i);      
end

numRows=length(Mat(:,1));
currRow=1;
% if Mat(1,1) == 1
%    RunLen=find(Mat(1,:),1,'last');
%    LB=1;
%    RB=align2(Mat(1,RunLen)+1);
%    Span=RB-LB;
%    Spacing=Span/(RunLen+1);
%    align2(Mat(1,1:RunLen))=align2(LB)+(1:RunLen)*Spacing;
%    currRow=currRow+1;
% else
%    while currRow <= numRows
%       
%    end
% end
a1=align1;
a2=align2;
while currRow <= numRows
   RunLen=find(Mat(currRow,:),1,'last');  
   if currRow==1 && Mat(1,1)==1
      LB=1;
   else
      LB=align2(Mat(currRow,1)-1);
   end
   if currRow==numRows && Mat(currRow,RunLen)==length(align2)
      RB=length(Indices2);
   else
      
      RB=align2(Mat(currRow,RunLen)+1);
   end
   Span=RB-LB;
   if Span == 0
      align2(Mat(currRow,1:RunLen))=LB;
   else
      Spacing=Span/(RunLen+1);
      align2(Mat(currRow,1:RunLen))=LB+(1:RunLen)*Spacing; 
   end
   currRow=currRow+1;
end
align2=round(align2);
for p=GappedAs
  tmp = 0;
  ct = 0;
  while tmp == 0
    if tmp == 0
      tmp = align2(max(p-ct,1));
    end
    if tmp == 0
      tmp = align2(min(p+ct,length(align2)));
    end
    ct=ct+1;
  end
  align2(p)=tmp;
end
