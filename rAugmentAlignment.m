% rAugmentAlignment(align1,align2) fills in a near match for non-aligned nucleotides
% For each nt in A aligned with a gap, the nearest aligned nt in B is found
% to use as the center of the band for that nt

function [align1, align2] = rAugmentAlignment(Indices1,align1,align2)

for i = 1:length(Indices1),
  if ~any(align1==i)
    align1 = [align1 i]; %#ok<AGROW>
    align2 = [align2 0]; %#ok<AGROW>
  end
end
[S1 IX] = sort(align1);
align1 = S1; %#ok<NASGU>
align2 = align2(IX);

GappedAs=find(align2==0);
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
