clear;
addpath([pwd filesep 'R3DAlign']);
Band{1}=20;
Band{2}=[40 10];
Band{3}=[200 70 20];
Neigh{1}=5;
Neigh{2}=[1 3];
Neigh{3}=[1 3 9];

%%%%%%%%%%%5S alignments%%%%%%%%%%%%%%%%%%
% clearvars -except Band Neigh
clearex('Band','Neigh');
Molecule{1} = '2QBG';   %E Coli
Molecule{2} = '2ZJR';   %Deinococcus Radiodurans 2ZJR
Chain{1}='A';   %5S
Chain{2}='Y';   %5S

clear NTList1 NTList2;
NTList1{1}='all';
NTList2{1}='all';
Query.Type='local';
Query.LoadFinal=0;
d=.4;
  for m1=1:length(Molecule)-1
    for m2=m1+1:length(Molecule) 
%       clearvars -except Band Neigh Molecule Chain NTList1 NTList2 Query d
%       m1 m2;
      clearex ('Band','Neigh','Molecule','Chain','NTList1','NTList2','Query','d','m1','m2');
      [Molecule{m1} ' ' Molecule{m2}]
      for iter=1:3
        for i=1:length(Band{iter}(:,1))
          for j=1:length(Neigh{iter}(:,1))
            clear B P D;
            for k=1:iter
              B{k}=Band{iter}(i,k);
              P{k}=Neigh{iter}(j,k);
              D{k}=d;
              CM{k}='greedy';
            end
            clear Chain1 Chain2;
            Chain1{1}=Chain{m1};
            Chain2{1}=Chain{m2};
            R3DAlign(Molecule{m1},Chain1,NTList1,Molecule{m2},Chain2,NTList2,D,P,B,CM,Query);
          end
        end
      end
    end      
  end
