% clear all;
Band{1}=[10;30;60;100];
Band{2}=[200 60;
         200 40;
         200 20;
         100 60;
%          100 40;
         100 20;
         100 10;
         60 40;
%          60 20;
         60 10;
         40 20;
         40 10];
Band{3}=[200 100 20;
%          200 100 50;
%          200 60 20;
         200 40 10;
         100 60 20;
         100 40 10;
%          60 40 20;
         60 40 10;
%          60 30 5;
         40 20 5];
% Band{4}=[200 100 50 10;
%          100 80 40 15;
%          100 60 25 10;
%          60 40 20 5];
Neigh{1}=[2;3;5;7;9];
Neigh{2}=[1 5;
           1 7;
           2 4;
%            2 6;
           2 8;
           2 10;
           3 5 ;
%            3 7;
           3 9;
           5 7;
%            5 9;
           5 11];
Neigh{3}=[1 4 8;
          1 3 9;
          2 5 8;
%           2 4 6;
          2 3 4;
          3 5 9;
%           3 7 9;
          3 4 5;
          4 6 8];
% Neigh{4}=[1 2 3 4;
%            1 3 5 7;
%            1 2 4 7;
%            1 3 5 9];

% %%%%%%%%%16S alignments%%%%%%%%%%%%%%%%%%       
% % % % % Molecule{1} = '1J5E';
% % % % % Molecule{2} = '2AW7';
% % % % % Molecule{3} = '2XZM';
% % % % % Molecule{4} = '3O30';
% % % % % Molecule{5} = '3JYV';
% % % % % Chain{1}='A';
% % % % % Chain{2}='A';
% % % % % Chain{3}='A';
% % % % % Chain{4}='1';
% % % % % Chain{5}='A';         
% Molecule{1} = '1FJG';
% Molecule{2} = '2AW7';
% Molecule{3} = '2XZM';
% Molecule{4} = '3U5F';
% Molecule{5} = '3JYV';
% Chain{1}='A';
% Chain{2}='A';
% Chain{3}='A';
% Chain{4}='6';
% Chain{5}='A';
% NTList1{1}='all';
% NTList2{1}='all';
% Query.Type='local';
% 
% for d=[.5 .5]
% %   for m1=1:length(Molecule)-1
% for m1=2
%     for m2=m1+1:length(Molecule) 
%     clearvars -except Band Neigh Molecule Chain NTList1 NTList2 Query d m1 m2;
%       [Molecule{m1} ' ' Molecule{m2}]
%       for iter=1:3
%         for i=1:length(Band{iter})
%           for j=1:length(Neigh{iter})
%             clear B P D;
%             for k=1:iter
%               B{k}=Band{iter}(i,k);
%               P{k}=Neigh{iter}(j,k);
%               D{k}=d;
%               CM{k}='greedy';
%             end
%             clear Chain1 Chain2;
%             Chain1{1}=Chain{m1};
%             Chain2{1}=Chain{m2};
%             D
%             B
%             P
%             R3DAlign(Molecule{m1},Chain1,NTList1,Molecule{m2},Chain2,NTList2,D,P,B,CM,Query);
%           end
%         end
%       end
%     end      
%   end
% end
% %%%%%%%%%%%23S alignments%%%%%%%%%%%%%%%%%%
% clearvars -except Band Neigh
% % % % % % % Molecule{1} = '2QBG';   %E Coli, Chain B
% % % % % % % Molecule{2} = '2ZJR';   %Deinococcus Radiodurans 2ZJR, Chain X
% % % % % % % Molecule{3} = '3PYO';   %Thermus Thermophilus 3PYO Chain A
% % % % % % % Molecule{4} = '1S72';   %Haloarcula marismortui 1S72, Chain 0
% Molecule{1} = '2QBG';   %E Coli, Chain B
% Molecule{2} = '2ZJR';   %Deinococcus Radiodurans 2ZJR, Chain X
% Molecule{3} = '3V2F';   %Thermus Thermophilus 3PYO Chain A
% Molecule{4} = '1S72';   %Haloarcula marismortui 1S72, Chain 0
% Chain{1}='B';
% Chain{2}='X';
% Chain{3}='A';
% Chain{4}='0';
% clear NTList1 NTList2;
% NTList1{1}='all';
% NTList2{1}='all';
% Query.Type='local';
% for d=[.4 .5]
%   for m1=1:length(Molecule)-1
%     for m2=m1+1:length(Molecule) 
%       clearvars -except Band Neigh Molecule Chain NTList1 NTList2 Query d m1 m2;
%       [Molecule{m1} ' ' Molecule{m2}]
%       for iter=1:3
%         for i=1:length(Band{iter})
%           for j=1:length(Neigh{iter})
%             clear B P D;
%             for k=1:iter
%               B{k}=Band{iter}(i,k);
%               P{k}=Neigh{iter}(j,k);
%               D{k}=d;
%               CM{k}='greedy';
%             end
%             clear Chain1 Chain2;
%             Chain1{1}=Chain{m1};
%             Chain2{1}=Chain{m2};
%             D
%             B
%             P
%             R3DAlign(Molecule{m1},Chain1,NTList1,Molecule{m2},Chain2,NTList2,D,P,B,CM,Query);
%           end
%         end
%       end
%     end      
%   end
% end
% % 3O58 is handled here because it includes 2 chains
% Molecule{5} = '3O58';   %Saccharomyces Cerevisiae 4 Angstrom 3O58 Chains 3, 1
% clear NTList1 NTList2;
% clear Chain1 Chain2;
% Chain2{1}='3';Chain2{2}='1';
% NTList1{1}='all';
% NTList2{1}='all';
% NTList2{2}='all';
% for d=[.4 .5]
%   for m1=1:length(Molecule)-1
%     clearvars -except Band Neigh Molecule Chain NTList1 NTList2 Query d m1 m2;
%     [Molecule{m1} ' ' Molecule{5}]
%     for iter=1:3
%       for i=1:length(Band{iter})
%         for j=1:length(Neigh{iter})
%           clear B P D;
%           for k=1:iter
%             B{k}=Band{iter}(i,k);
%             P{k}=Neigh{iter}(j,k);
%             D{k}=d;
%             CM{k}='greedy';
%           end
%           clear Chain1;
%           Chain1{1}=Chain{m1};
%           D
%           P
%           B
%           R3DAlign(Molecule{m1},Chain1,NTList1,Molecule{5},Chain2,NTList2,D,P,B,CM,Query);
%         end
%       end
%     end
%   end      
% end

%%%%%%%%%%%5S alignments%%%%%%%%%%%%%%%%%%
clearvars -except Band Neigh
% % % % % % Molecule{1} = '2QBG';   %E Coli, Chain B
% % % % % % Molecule{2} = '2ZJR';   %Deinococcus Radiodurans 2ZJR, Chain X
% % % % % % Molecule{3} = '3PYO';   %Thermus Thermophilus 3PYO Chain A
% % % % % % Molecule{4} = '1S72';   %Haloarcula marismortui 1S72, Chain 0
Molecule{1} = '2QBG';   %E Coli
Molecule{2} = '2ZJR';   %Deinococcus Radiodurans 2ZJR
Molecule{3} = '3V2F';   %Thermus Thermophilus 3PYO 
Molecule{4} = '1S72';   %Haloarcula marismortui 1S72
Molecule{5} = '3O58';   %Saccharomyces Cerevisiae 4 Angstrom
Chain{1}='A';   %5S
Chain{2}='Y';   %5S
Chain{3}='B';   %5S
Chain{4}='9';   %5S
Chain{5}='2';   %5S
clear NTList1 NTList2;
NTList1{1}='all';
NTList2{1}='all';
Query.Type='local';
for d=[.4 .5]
  for m1=1:length(Molecule)-1
    for m2=m1+1:length(Molecule) 
      clearvars -except Band Neigh Molecule Chain NTList1 NTList2 Query d m1 m2;
      [Molecule{m1} ' ' Molecule{m2}]
      for iter=1:3
        for i=1:length(Band{iter})
          for j=1:length(Neigh{iter})
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
            D
            B
            P
            R3DAlign(Molecule{m1},Chain1,NTList1,Molecule{m2},Chain2,NTList2,D,P,B,CM,Query);
          end
        end
      end
    end      
  end
end
