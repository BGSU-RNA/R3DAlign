%rFindAlignmentDiscrepancies assigns a discrepancy for each nucleotide in
%an alignment according to one of two possible methods:
%  1) For each nt in structure A, it and the four nearest nt's in A are
%     superimposed onto the corresponding nt's in B
%  2) For each nt in A, the 10 neighborhoods containing it and 3 other
%     nucleotides with the smallest diameter are superimposed onto the
%     corresponding nt's in B.  The median of the 10 discrepancies is used.

function [D,GoodDiscreps] = rFindAlignmentDiscrepancies(File1,Indices1,File2,Indices2,Method)

c1 = cat(1,File1.NT(Indices1).Center); % nucleotide centers
D1 = full(zMutualDistance(c1,Inf));

GoodDiscreps=[];

if isequal(Method,'nearest4')
   Discreps4=zeros(1,length(Indices1));
   for i=1:length(Indices1)
      [a b]=sort(D1(i,:));
      Discreps4(i) = xDiscrepancy(File1,Indices1(b(1:5)),File2,Indices2(b(1:5)));
      if Discreps4(i)<.2
          ItoAdd=[Indices1(b(1:5))' Indices2(b(1:5))'];
          GoodDiscreps=[GoodDiscreps; ItoAdd]; %#ok<AGROW>
      end
   end
   D=Discreps4';
elseif isequal(Method,'AllQuads')
   A = triu(D1);
   [rA cA] = find(A<10 & A>0);
   [S1 IX] = sort(rA);
   rA = S1;
   cA = cA(IX);
   ATrips = DoublesToTriples([rA cA]);
   AQuads = TriplesToQuads(ATrips);
   AQuads = rRefineQuads(AQuads,A,8);
   nA=length(AQuads(:,1));
   Discrepancies=cell(length(Indices1),1);

   for i = 1:nA
      Disc = xDiscrepancy(File1,Indices1(AQuads(i,:)),File2,Indices2(AQuads(i,:)));
      Discrepancies{AQuads(i,1),1}=[Discrepancies{AQuads(i,1),1} Disc];
      Discrepancies{AQuads(i,2),1}=[Discrepancies{AQuads(i,2),1} Disc];
      Discrepancies{AQuads(i,3),1}=[Discrepancies{AQuads(i,3),1} Disc];
      Discrepancies{AQuads(i,4),1}=[Discrepancies{AQuads(i,4),1} Disc];
   end

   DiscrepsQ=zeros(1,length(Indices1));
   for i = 1:length(Indices1)
      DiscrepsQ(i)=median(Discrepancies{i});
   end
   D=DiscrepsQ';
else
   fprintf('ERROR in rFindAlignmentDiscrepancies - invalid Method entered\n')
end