%This function takes a matrix (either 3 by N or N by 3) of satisfactory
%triples and constructs a list of triples (TripCt by 3 matrix) so that all
%triples are in triples list

function[QuadsList] = TriplesToQuads(TriplesList)

TL = TriplesList;
QuadsList=[];
if length(TL(:,1))==3
    TL=TL';
end

for m = transpose(unique(TL(:,1)))
    I3 = find(TL(:,1)==m);
    %now just find a triple in that subset, add m to it, and there's a quad
    mTrips=DoublesToTriples([TL(I3,2) TL(I3,3)]);
    if ~isempty(mTrips)
       newQuads = horzcat(m*ones(length(mTrips(:,1)),1),mTrips);
       QuadsList = vertcat(QuadsList,newQuads);
    end
end
