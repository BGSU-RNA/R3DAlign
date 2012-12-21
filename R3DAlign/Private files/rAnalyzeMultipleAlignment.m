% Evaluate the number of near basepairs aligned with  true basepairs, near with near, 
% true with true, etc.

function [BPFamily] = rAnalyzeMultipleAlignment(File1,File2,File3,SpreadsheetName)
[a, b]=xlsread(SpreadsheetName); %#ok<ASGLU>
for k = [2 5 8]
   for i=1:length(b(:,k))
      b(i,k)=lower(strtrim(b(i,k)));
      if isequal(b{i,k},'tsh')
         b{i,k}='ths';
      elseif isequal(b{i,k},'csh')
         b{i,k}='chs';
      elseif isequal(b{i,k},'tsw')
         b{i,k}='tws';
      elseif isequal(b{i,k},'csw')
         b{i,k}='cws';
      elseif isequal(b{i,k},'thw')
         b{i,k}='twh';
      elseif isequal(b{i,k},'chw')
         b{i,k}='cwh';
      elseif isequal(b{i,k},'ntsh')
         b{i,k}='nths';
      elseif isequal(b{i,k},'ncsh')
         b{i,k}='nchs';
      elseif isequal(b{i,k},'ntsw')
         b{i,k}='ntws';
      elseif isequal(b{i,k},'ncsw')
         b{i,k}='ncws';
      elseif isequal(b{i,k},'nthw')
         b{i,k}='ntwh';
      elseif isequal(b{i,k},'nchw')
         b{i,k}='ncwh';
      elseif ~isempty(b{i,k}) && isequal(b{i,k}(1),'s')
         b{i,k}='stacking';
      end
   end
end
Names={'cww','tww','css','tss','chh','thh','chs','ths','cws','tws','cwh','twh'};
for i=1:12
   BPFamily(i).Name=Names{i}; %#ok<*AGROW>
   BPFamily(i).numAligned=0;
   BPFamily(i).numAlignedConserved=0;
   BPFamily(i).freqsAll3=zeros(17,1);
   BPFamily(i).freqs12=zeros(17,17);
   BPFamily(i).freqs13=zeros(17,17);
   BPFamily(i).freqs23=zeros(17,17);
   BPFamily(i).nearfreqs12=zeros(17,17);
   BPFamily(i).nearfreqs13=zeros(17,17);
   BPFamily(i).nearfreqs23=zeros(17,17);
   BPFamily(i).relfreqs=zeros(17,17);
   BPFamily(i).nearnumAligned=0;
   BPFamily(i).nearfreqsAll3=zeros(17,1);
   BPFamily(i).nearrelfreqs=zeros(17,17);
   BPFamily(i).nearnumAlignedConserved=0;
   BPFamily(i).numAlignedConserved12not3=0;
   BPFamily(i).nearnumAlignedConserved12not3=0;
   BPFamily(i).numAlignedConserved13not2=0;
   BPFamily(i).nearnumAlignedConserved13not2=0;
   BPFamily(i).numAlignedConserved23not1=0;
   BPFamily(i).nearnumAlignedConserved23not1=0;
   BPFamily(i).numAlignedConservedNone=0;
   BPFamily(i).nearnumAlignedConservedNone=0;
end

for i=2:length(b(:,2))
   a1=b{i,1}(1);
   b1=b{i,3}(1);
   a2=b{i,4}(1);
   b2=b{i,6}(1);
   a3=b{i,7}(1);
   b3=b{i,9}(1);
   a1=char2num(a1);
   b1=char2num(b1);
   a2=char2num(a2);
   b2=char2num(b2);
   a3=char2num(a3);
   b3=char2num(b3);
   if a1~=0 && b1~=0 && zIndexLookup(File1,b{i,1}(2:end),'X') < zIndexLookup(File1,b{i,3}(2:end),'X')
      for j=1:12
         if isequal(strtrim(b{i,2}),Names{j})
            if isequal(b{i,5},Names{j})
               if isequal(b{i,8},Names{j})
                  BPFamily(j).numAligned=BPFamily(j).numAligned+1;
                  if a1==a2 && a1==a3 && b1==b2 && b1==b3
                     BPFamily(j).numAlignedConserved=BPFamily(j).numAlignedConserved+1;
                     BPFamily(j).freqsAll3((a1-1)*4+b1)=BPFamily(j).freqsAll3((a1-1)*4+b1)+1;
                     BPFamily(j).freqsAll3(17)=BPFamily(j).freqsAll3(17)+1;
                  elseif a1==a2 && b1==b2
                     BPFamily(j).numAlignedConserved12not3=BPFamily(j).numAlignedConserved12not3+1;
                     BPFamily(j).freqs12((a1-1)*4+b1,(a3-1)*4+b3)=BPFamily(j).freqs12((a1-1)*4+b1,(a3-1)*4+b3)+1;
                     BPFamily(j).freqs12((a1-1)*4+b1,17)=BPFamily(j).freqs12((a1-1)*4+b1,17)+1;
                     BPFamily(j).freqs12(17,(a3-1)*4+b3)=BPFamily(j).freqs12(17,(a3-1)*4+b3)+1;
                     BPFamily(j).freqs12(17,17)=BPFamily(j).freqs12(17,17)+1;
                  elseif a1==a3 && b1==b3
                     BPFamily(j).numAlignedConserved13not2=BPFamily(j).numAlignedConserved13not2+1;
                     BPFamily(j).freqs13((a1-1)*4+b1,(a2-1)*4+b2)=BPFamily(j).freqs13((a1-1)*4+b1,(a2-1)*4+b2)+1;
                     BPFamily(j).freqs13((a1-1)*4+b1,17)=BPFamily(j).freqs13((a1-1)*4+b1,17)+1;
                     BPFamily(j).freqs13(17,(a2-1)*4+b2)=BPFamily(j).freqs13(17,(a2-1)*4+b2)+1;
                     BPFamily(j).freqs13(17,17)=BPFamily(j).freqs13(17,17)+1;
                  elseif a2==a3 && b2==b3
                     BPFamily(j).numAlignedConserved23not1=BPFamily(j).numAlignedConserved23not1+1;
                     BPFamily(j).freqs23((a2-1)*4+b2,(a1-1)*4+b1)=BPFamily(j).freqs23((a2-1)*4+b2,(a1-1)*4+b1)+1;
                     BPFamily(j).freqs23((a2-1)*4+b2,17)=BPFamily(j).freqs23((a2-1)*4+b2,17)+1;
                     BPFamily(j).freqs23(17,(a1-1)*4+b1)=BPFamily(j).freqs23(17,(a1-1)*4+b1)+1;
                     BPFamily(j).freqs23(17,17)=BPFamily(j).freqs23(17,17)+1;
                  elseif (a1~=a2 || b1~=b2) && (a1~=a3 || b1~=b3) && (a2~=a3 || b2~=b3)
                     BPFamily(j).numAlignedConservedNone=BPFamily(j).numAlignedConservedNone+1;
                  end
               end
            end
         end
         if isequal(strtrim(b{i,2}),Names{j}) || isequal(strtrim(b{i,2}),['n' Names{j}]) 
            if isequal(b{i,5},Names{j}) || isequal(strtrim(b{i,5}),['n' Names{j}])
               if isequal(b{i,8},Names{j}) || isequal(strtrim(b{i,8}),['n' Names{j}])
                  BPFamily(j).nearnumAligned=BPFamily(j).nearnumAligned+1;
                  if a1==a2 && a1==a3 && b1==b2 && b1==b3
                     BPFamily(j).nearnumAlignedConserved=BPFamily(j).nearnumAlignedConserved+1;
                     BPFamily(j).nearfreqsAll3((a1-1)*4+b1)=BPFamily(j).nearfreqsAll3((a1-1)*4+b1)+1;
                     BPFamily(j).nearfreqsAll3(17)=BPFamily(j).nearfreqsAll3(17)+1;
                  elseif a1==a2 && b1==b2
                     BPFamily(j).nearnumAlignedConserved12not3=BPFamily(j).nearnumAlignedConserved12not3+1;
                     BPFamily(j).nearfreqs12((a1-1)*4+b1,(a3-1)*4+b3)=BPFamily(j).nearfreqs12((a1-1)*4+b1,(a3-1)*4+b3)+1;
                     BPFamily(j).nearfreqs12((a1-1)*4+b1,17)=BPFamily(j).nearfreqs12((a1-1)*4+b1,17)+1;
                     BPFamily(j).nearfreqs12(17,(a3-1)*4+b3)=BPFamily(j).nearfreqs12(17,(a3-1)*4+b3)+1;
                     BPFamily(j).nearfreqs12(17,17)=BPFamily(j).nearfreqs12(17,17)+1;
                  elseif a1==a3 && b1==b3
                     BPFamily(j).nearnumAlignedConserved13not2=BPFamily(j).nearnumAlignedConserved13not2+1;
                     BPFamily(j).nearfreqs13((a1-1)*4+b1,(a2-1)*4+b2)=BPFamily(j).nearfreqs13((a1-1)*4+b1,(a2-1)*4+b2)+1;
                     BPFamily(j).nearfreqs13((a1-1)*4+b1,17)=BPFamily(j).nearfreqs13((a1-1)*4+b1,17)+1;
                     BPFamily(j).nearfreqs13(17,(a2-1)*4+b2)=BPFamily(j).nearfreqs13(17,(a2-1)*4+b2)+1;
                     BPFamily(j).nearfreqs13(17,17)=BPFamily(j).nearfreqs13(17,17)+1;
                  elseif a2==a3 && b2==b3
                     BPFamily(j).nearnumAlignedConserved23not1=BPFamily(j).nearnumAlignedConserved23not1+1;
                     BPFamily(j).nearfreqs23((a2-1)*4+b2,(a1-1)*4+b1)=BPFamily(j).nearfreqs23((a2-1)*4+b2,(a1-1)*4+b1)+1;
                     BPFamily(j).nearfreqs23((a2-1)*4+b2,17)=BPFamily(j).nearfreqs23((a2-1)*4+b2,17)+1;
                     BPFamily(j).nearfreqs23(17,(a1-1)*4+b1)=BPFamily(j).nearfreqs23(17,(a1-1)*4+b1)+1;
                     BPFamily(j).nearfreqs23(17,17)=BPFamily(j).nearfreqs23(17,17)+1;
                  elseif (a1~=a2 || b1~=b2) && (a1~=a3 || b1~=b3) && (a2~=a3 || b2~=b3)
                     BPFamily(j).nearnumAlignedConservedNone=BPFamily(j).nearnumAlignedConservedNone+1;
                  end
               end
            end
         end
      end
   end
end
% % % for k=1:12
% % %    for i=1:16
% % %       for j=1:16
% % %        BPFamily(k).relfreqs(i,j)=round2(BPFamily(k).freqs(i,j)/BPFamily(k).freqs(i,17),.01);
% % %        BPFamily(k).nearrelfreqs(i,j)=round2(BPFamily(k).nearfreqs(i,j)/BPFamily(k).nearfreqs(i,17),.01);
% % %       end
% % %    end
% % %    BPFamily(k).relfreqs(:,17)=BPFamily(k).freqs(:,17);
% % %    BPFamily(k).relfreqs(17,:)=BPFamily(k).freqs(17,:);
% % %    BPFamily(k).nearrelfreqs(:,17)=BPFamily(k).nearfreqs(:,17);
% % %    BPFamily(k).nearrelfreqs(17,:)=BPFamily(k).nearfreqs(17,:);
% % % end       

% Label{1,2}='AA';Label{1,3}='AC';Label{1,4}='AG';Label{1,5}='AU';
% Label{1,6}='CA';Label{1,7}='CC';Label{1,8}='CG';Label{1,9}='CU';
% Label{1,10}='GA';Label{1,11}='GC';Label{1,12}='GG';Label{1,13}='GU';
% Label{1,14}='UA';Label{1,15}='UC';Label{1,16}='UG';Label{1,17}='UU';
% Label{2,1}='AA';Label{3,1}='AC';Label{4,1}='AG';Label{5,1}='AU';
% Label{6,1}='CA';Label{7,1}='CC';Label{8,1}='CG';Label{9,1}='CU';
% Label{10,1}='GA';Label{11,1}='GC';Label{12,1}='GG';Label{13,1}='GU';
% Label{14,1}='UA';Label{15,1}='UC';Label{16,1}='UG';Label{17,1}='UU';
% Label{17,17}='';

Label{1}='AA';Label{2}='AC';Label{3}='AG';Label{4}='AU';
Label{5}='CA';Label{6}='CC';Label{7}='CG';Label{8}='CU';
Label{9}='GA';Label{10}='GC';Label{11}='GG';Label{12}='GU';
Label{13}='UA';Label{14}='UC';Label{15}='UG';Label{16}='UU';
ExcelName=[File1.Filename '_' File2.Filename '_' File3.Filename '_alignment_data.xlsx'];
P1={};
P2={};
for i=1:12
   if mod(i,2)==1
      P2{(i-1)/2*19+1,1}=Names{i};
      P2((i-1)/2*19+1,2:17)=Label(1:16);
      P2((i-1)/2*19+2:(i-1)/2*19+17,1)=Label(1:16);
   else
      P2{(i-2)/2*19+1,20}=Names{i};
      P2((i-2)/2*19+1,21:36)=Label(1:16);
      P2((i-2)/2*19+2:(i-2)/2*19+17,20)=Label(1:16);
   end
end
P1(2:17,1)=Label(1:16);
P1{18,1}='TOTAL';
P3=P1;
BPFamily=Reformat(BPFamily)
for i=1:12
   P1{1,i+1}=Names{i};
   P3{1,i+1}=Names{i};
   P1(2:18,i+1)=BPFamily(i).freqsAll3;
   P3(2:18,i+1)=BPFamily(i).nearfreqsAll3;
end

% P4=P2;
P5=P2;
P6=P2;
P7=P2;
P8=P2;
P9=P2;
P10=P2;
for i=1:12
   if mod(i,2)==1
%       P1((i-1)/2*19+2,2:18)=num2cell(BPFamily(i).freqsAll3);
%       P2((i-1)/2*19+2:(i-1)/2*19+18,2:18)=num2cell(BPFamily(i).relfreqs);
%       P3((i-1)/2*19+2,2:18)=num2cell(BPFamily(i).nearfreqsAll3);
%       P4((i-1)/2*19+2:(i-1)/2*19+18,2:18)=num2cell(BPFamily(i).nearrelfreqs);
      P5((i-1)/2*19+2:(i-1)/2*19+18,2:18)=BPFamily(i).freqs12;
      P6((i-1)/2*19+2:(i-1)/2*19+18,2:18)=BPFamily(i).freqs13;
      P7((i-1)/2*19+2:(i-1)/2*19+18,2:18)=BPFamily(i).freqs23;
      P8((i-1)/2*19+2:(i-1)/2*19+18,2:18)=BPFamily(i).nearfreqs12;
      P9((i-1)/2*19+2:(i-1)/2*19+18,2:18)=BPFamily(i).nearfreqs13;
      P10((i-1)/2*19+2:(i-1)/2*19+18,2:18)=BPFamily(i).nearfreqs23;
   else
%       P1((i-2)/2*19+2,21:37)=num2cell(BPFamily(i).freqsAll3);
%       P2((i-2)/2*19+2:(i-2)/2*19+18,21:37)=num2cell(BPFamily(i).relfreqs);
%       P3((i-2)/2*19+2,21:37)=num2cell(BPFamily(i).nearfreqsAll3);
%       P4((i-2)/2*19+2:(i-2)/2*19+18,21:37)=num2cell(BPFamily(i).nearrelfreqs);
      P5((i-2)/2*19+2:(i-2)/2*19+18,21:37)=BPFamily(i).freqs12;
      P6((i-2)/2*19+2:(i-2)/2*19+18,21:37)=BPFamily(i).freqs13;
      P7((i-2)/2*19+2:(i-2)/2*19+18,21:37)=BPFamily(i).freqs23;
      P8((i-2)/2*19+2:(i-2)/2*19+18,21:37)=BPFamily(i).nearfreqs12;
      P9((i-2)/2*19+2:(i-2)/2*19+18,21:37)=BPFamily(i).nearfreqs13;
      P10((i-2)/2*19+2:(i-2)/2*19+18,21:37)=BPFamily(i).nearfreqs23;
   end
end

xlswrite(ExcelName,P1,'All3','A1');
% xlswrite(ExcelName,P2,'BPCovRelFreqs','A1');
xlswrite(ExcelName,P3,'All3Near','A1');
% xlswrite(ExcelName,P4,'BPCovRelFreqsNear','A1');
xlswrite(ExcelName,P5,'12','A1');
xlswrite(ExcelName,P6,'13','A1');
xlswrite(ExcelName,P7,'23','A1');
xlswrite(ExcelName,P8,'12Near','A1');
xlswrite(ExcelName,P9,'13Near','A1');
xlswrite(ExcelName,P10,'23Near','A1');
% % % for i=1:12
% % %    Label{1,1}=Names{i};
% % %    if mod(i,2)==1
% % %       xlswrite(ExcelName,Label,'BPCovFreqs',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
% % %       xlswrite(ExcelName,BPFamily(i).freqs,'BPCovFreqs',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
% % %       xlswrite(ExcelName,Label,'BPCovRelFreqs',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
% % %       xlswrite(ExcelName,BPFamily(i).relfreqs,'BPCovRelFreqs',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
% % %       xlswrite(ExcelName,Label,'BPCovFreqsNear',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
% % %       xlswrite(ExcelName,BPFamily(i).nearfreqs,'BPCovFreqsNear',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
% % %       xlswrite(ExcelName,Label,'BPCovRelFreqsNear',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
% % %       xlswrite(ExcelName,BPFamily(i).nearrelfreqs,'BPCovRelFreqsNear',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
% % %    else
% % %       xlswrite(ExcelName,Label,'BPCovFreqs',['T' num2str((i-2)/2*19+1) ':AJ' num2str((i-2)/2*19+17)]);
% % %       xlswrite(ExcelName,BPFamily(i).freqs,'BPCovFreqs',['U' num2str((i-2)/2*19+2) ':AK' num2str((i-2)/2*19+18)]);
% % %       xlswrite(ExcelName,Label,'BPCovRelFreqs',['T' num2str((i-2)/2*19+1) ':AJ' num2str((i-2)/2*19+17)]);
% % %       xlswrite(ExcelName,BPFamily(i).relfreqs,'BPCovRelFreqs',['U' num2str((i-2)/2*19+2) ':AK' num2str((i-2)/2*19+18)]); 
% % %       xlswrite(ExcelName,Label,'BPCovFreqsNear',['T' num2str((i-2)/2*19+1) ':AJ' num2str((i-2)/2*19+17)]);
% % %       xlswrite(ExcelName,BPFamily(i).nearfreqs,'BPCovFreqsNear',['U' num2str((i-2)/2*19+2) ':AK' num2str((i-2)/2*19+18)]);
% % %       xlswrite(ExcelName,Label,'BPCovRelFreqsNear',['T' num2str((i-2)/2*19+1) ':AJ' num2str((i-2)/2*19+17)]);
% % %       xlswrite(ExcelName,BPFamily(i).nearrelfreqs,'BPCovRelFreqsNear',['U' num2str((i-2)/2*19+2) ':AK' num2str((i-2)/2*19+18)]); 
% % %    end
% % % end

S{2,1}=['Basepairs conserved in ' File1.Filename ' and ' File2.Filename ' and ' File3.Filename];
S{3,1}='Basepairs and base combinations conserved in all 3';
S{4,1}=['Basepairs conserved in all 3 and base combinations conserved in ' File1.Filename ' and ' File2.Filename ' but not ' File3.Filename];
S{5,1}=['Basepairs conserved in all 3 and base combinations conserved in ' File1.Filename ' and ' File3.Filename ' but not ' File2.Filename];
S{6,1}=['Basepairs conserved in all 3 and base combinations conserved in ' File2.Filename ' and ' File3.Filename ' but not ' File1.Filename];
S{7,1}='Basepairs conserved in all 3 and base combinations conserved in none';
% S{8,1}='Percentage of basepairs with conserved bases';
S{9,1}='Basepairs conserved in all 3 (counting near-bps)';
S{10,1}=['Basepairs and bases conserved in ' File1.Filename ' and ' File2.Filename ' and ' File3.Filename ' (counting near-bps)'];
S{11,1}=['Basepairs conserved in all 3 and base combinations conserved in ' File1.Filename ' and ' File2.Filename ' but not ' File3.Filename ' (counting near-bps)'];
S{12,1}=['Basepairs conserved in all 3 and base combinations conserved in ' File1.Filename ' and ' File3.Filename ' but not ' File2.Filename ' (counting near-bps)'];
S{13,1}=['Basepairs conserved in all 3 and base combinations conserved in ' File2.Filename ' and ' File3.Filename ' but not ' File1.Filename ' (counting near-bps)'];
S{14,1}='Basepairs conserved in all 3 and base combinations conserved in none (counting near-bps)';
% S{15,1}='Percentage of basepairs with conserved bases';
% S(1,2:13)=Names;
for i=1:12
   S{1,2*i}=Names{i};
   S{2,2*i}=BPFamily(i).numAligned;
   S{3,2*i}=BPFamily(i).numAlignedConserved;
   S{3,2*i+1}=round2(BPFamily(i).numAlignedConserved/BPFamily(i).numAligned,.01)*100;
   S{3,2*i+1}=[num2str(S{3,2*i+1}) '%'];
   S{4,2*i}=BPFamily(i).numAlignedConserved12not3;
   S{4,2*i+1}=round2(BPFamily(i).numAlignedConserved12not3/BPFamily(i).numAligned,.01)*100;
   S{4,2*i+1}=[num2str(S{4,2*i+1}) '%'];
   S{5,2*i}=BPFamily(i).numAlignedConserved13not2;
   S{5,2*i+1}=round2(BPFamily(i).numAlignedConserved13not2/BPFamily(i).numAligned,.01)*100;
   S{5,2*i+1}=[num2str(S{5,2*i+1}) '%'];
   S{6,2*i}=BPFamily(i).numAlignedConserved23not1;
   S{6,2*i+1}=round2(BPFamily(i).numAlignedConserved23not1/BPFamily(i).numAligned,.01)*100;
   S{6,2*i+1}=[num2str(S{6,2*i+1}) '%'];
   S{7,2*i}=BPFamily(i).numAlignedConservedNone;
   S{7,2*i+1}=round2(BPFamily(i).numAlignedConservedNone/BPFamily(i).numAligned,.01)*100;
   S{7,2*i+1}=[num2str(S{7,2*i+1}) '%'];
   
   S{9,2*i}=BPFamily(i).nearnumAligned;
   S{10,2*i}=BPFamily(i).nearnumAlignedConserved;
   S{10,2*i+1}=round2(BPFamily(i).nearnumAlignedConserved/BPFamily(i).nearnumAligned,.01)*100;
   S{10,2*i+1}=[num2str(S{10,2*i+1}) '%'];
   S{11,2*i}=BPFamily(i).nearnumAlignedConserved12not3;
   S{11,2*i+1}=round2(BPFamily(i).nearnumAlignedConserved12not3/BPFamily(i).nearnumAligned,.01)*100;
   S{11,2*i+1}=[num2str(S{11,2*i+1}) '%'];
   S{12,2*i}=BPFamily(i).nearnumAlignedConserved13not2;
   S{12,2*i+1}=round2(BPFamily(i).nearnumAlignedConserved13not2/BPFamily(i).nearnumAligned,.01)*100;
   S{12,2*i+1}=[num2str(S{12,2*i+1}) '%'];
   S{13,2*i}=BPFamily(i).nearnumAlignedConserved23not1;
   S{13,2*i+1}=round2(BPFamily(i).nearnumAlignedConserved23not1/BPFamily(i).nearnumAligned,.01)*100;
   S{13,2*i+1}=[num2str(S{13,2*i+1}) '%'];
   S{14,2*i}=BPFamily(i).nearnumAlignedConservedNone;
   S{14,2*i+1}=round2(BPFamily(i).nearnumAlignedConservedNone/BPFamily(i).nearnumAligned,.01)*100;
   S{14,2*i+1}=[num2str(S{14,2*i+1}) '%'];
  
end
xlswrite(ExcelName,S,'Summary','A1');

return;
end

function [num]=char2num(letter)
  letter = lower(letter);
  if isequal(letter,'a')
     num=1;
  elseif isequal(letter,'c')
     num=2;
  elseif isequal(letter,'g')
     num=3;
  elseif isequal(letter,'u')
     num=4;
  else
     num=0;
  end
end

function [charcell]=nummat2charcell(A)
   [m n]=size(A);
   charcell=cell(m,n);
   for i=1:m
      for j=1:n
         charcell{i,j}=num2str(A(i,j));
      end
   end
end

function [BPF]=Reformat(BPF)
   for i=1:12
      BPF(i).freqsAll3=nummat2charcell(BPF(i).freqsAll3);
      BPF(i).freqs12=nummat2charcell(BPF(i).freqs12);
      BPF(i).freqs13=nummat2charcell(BPF(i).freqs13);
      BPF(i).freqs23=nummat2charcell(BPF(i).freqs23);
      BPF(i).nearfreqsAll3=nummat2charcell(BPF(i).nearfreqsAll3);
      BPF(i).nearfreqs12=nummat2charcell(BPF(i).nearfreqs12);
      BPF(i).nearfreqs13=nummat2charcell(BPF(i).nearfreqs13);
      BPF(i).nearfreqs23=nummat2charcell(BPF(i).nearfreqs23);
      for j=1:16
            if isequal(BPF(i).freqsAll3{j,1},'0')
               BPF(i).freqsAll3{j,1}='';
            end
            if isequal(BPF(i).nearfreqsAll3{j,1},'0')
               BPF(i).nearfreqsAll3{j,1}='';
            end
      end
      for j=1:16
         for k=1:16
            if isequal(BPF(i).freqs12{j,k},'0')
               BPF(i).freqs12{j,k}='';
            end
            if isequal(BPF(i).nearfreqs12{j,k},'0')
               BPF(i).nearfreqs12{j,k}='';
            end
            if isequal(BPF(i).freqs13{j,k},'0')
               BPF(i).freqs13{j,k}='';
            end
            if isequal(BPF(i).nearfreqs13{j,k},'0')
               BPF(i).nearfreqs13{j,k}='';
            end
            if isequal(BPF(i).freqs23{j,k},'0')
               BPF(i).freqs23{j,k}='';
            end
            if isequal(BPF(i).nearfreqs23{j,k},'0')
               BPF(i).nearfreqs23{j,k}='';
            end
         end
      end
   end
end
% bothblank+match+nearmatch+notmatch
% samebp+diffbp+inAnotB+AwithNear+AwithWrongNear+AwithNoNT+AwithOneNT+inBnotA...
% +BwithNear+BwithWrongNear+BwithNoNT+BwithOneNT+SingleAwithNot+SingleBwithNot+SinglesAligned;
