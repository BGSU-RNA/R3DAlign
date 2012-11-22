function [newQL] = rGetQuads(QL,DM,refineNum)     %QuadsList and Distance Matrix

N=length(DM);
newQL=[];
numNeigh=zeros(N,2);
numNeigh(1:N,1)=1:N;

for i=1:length(QL(:,1))
   numNeigh(QL(i,1),2)=numNeigh(QL(i,1),2)+1;
   numNeigh(QL(i,2),2)=numNeigh(QL(i,2),2)+1;
   numNeigh(QL(i,3),2)=numNeigh(QL(i,3),2)+1;
   numNeigh(QL(i,4),2)=numNeigh(QL(i,4),2)+1;
end

low=find(numNeigh(:,2)<refineNum);

if ~isequal(DM, DM')
   DM=DM+DM';
end

NUM=5;
while nchoosek(NUM,3) < refineNum
   NUM = NUM + 1;
end

for k = low'
    [y z]=sort(DM(k,:));
    T=z(2:NUM+1);
    for a = 1:(length(T)-2)
       for b = (a+1):(length(T)-1)
           for c = (b+1):length(T)
              Cand = sort([k T(a) T(b) T(c)]);
              [r2] = find(QL(:,1)==Cand(1));
              found=false;
              for m=r2'
                 if isequal(Cand,QL(m,:))
                    found=true;
                    break;
                 end
              end
              if found==false
                 numNeigh(Cand(1),2)=numNeigh(Cand(1),2)+1;
                 numNeigh(Cand(2),2)=numNeigh(Cand(2),2)+1;
                 numNeigh(Cand(3),2)=numNeigh(Cand(3),2)+1;
                 numNeigh(Cand(4),2)=numNeigh(Cand(4),2)+1;
                 QL=vertcat(QL, Cand);
              end
              
           end
       end
    end
end

QDiam=zeros(length(QL(:,1)),1);
for i=1:length(QL(:,1)) 
  QDiam(i) = max(max(DM(QL(i,:),QL(i,:))));
end

for i=1:N
   if numNeigh(i,2)>0
      [R C] = find(QL==i);
      clear C;
      T=QDiam(R);
      [V I] = sort(T);
%    eg. I(1:10) is the indices into R of the 10 neighborhoods with lowest
%    diameters
      I(1:min(numNeigh(i,2),refineNum));
      Cand = QL(R(I(1:min(numNeigh(i,2),refineNum))),:);    %the best quads containing NT i
      for j=1:length(Cand(:,1))
         if isempty(newQL)
            newQL=Cand(j,:);
         else
            [r1 c1] = find(newQL(:,1)==Cand(j,1)); %#ok<*NASGU>
            found=false;
            for k=r1'
               if isequal(Cand(j,:),newQL(k,:))
                  found=true;
                  break;
               end
            end
            if found==false
               newQL=vertcat(newQL,Cand(j,:));
            end
         end
      end
   end
end

newQL=sortrows(newQL);
