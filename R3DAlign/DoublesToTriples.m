%This function takes a matrix (either 2 by N or N by 2) of satisfactory
%pairs and constructs a list of triples (TripCt by 3 matrix) that satisfies all pairwise
%constraints

function[TriplesList] = DoublesToTriples(DoublesList)

if length(DoublesList(:,1)) > 1000
   TriplesList=zeros(100000,3);
else
   TriplesList=zeros(1000,3);
end

DL = DoublesList;

if length(DL(:,1))==2
   if length(DL(1,:))==2
      TriplesList=[];
      return;
   else
      DL=DL';
   end
end

N = length(DoublesList(:,1));
DM=sparse(max(max(DoublesList)),max(max(DoublesList)));
TripCt = 0;

for i = 1:N
   DM(DL(i,1),DL(i,2))=1;
end

for i = transpose(unique(DL(:,1)))
   I2 = find(DL(:,1)==i);
   for j = 1:length(I2)-1
       for k = j+1:length(I2)
           if DM(DL(I2(j),2),DL(I2(k),2)) == 1
               TripCt = TripCt + 1;
               TriplesList(TripCt,:) = [i DL(I2(j),2) DL(I2(k),2)];
           end
       end
   end
end
if TripCt==0
   TriplesList=[];
end
TriplesList=TriplesList(1:TripCt,:);