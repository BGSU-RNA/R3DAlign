% 
a=1;
Al(a).Name='Ryan812';
Al(a).SheetName='Ryan812';
a=a+1;
Al(a).Name='Ryan8132';
Al(a).SheetName='Ryan8132';
a=a+1;
Al(a).Name='SARA';
Al(a).SheetName='SARA';
a=a+1;
Al(a).Name='RyanSARA';
Al(a).SheetName='RyanSARA';
a=a+1;
Al(a).Name='RyanSARA2';
Al(a).SheetName='RyanSARA2';
% a=a+1;
% Al(a).Name='SARSA';
% Al(a).SheetName='SARSA';
% a=a+1;
% Al(a).Name='SARSA-Ryan';
% Al(a).SheetName='SARSA-Ryan';

T=cell(17,length(Al)+1);

for k = 1:length(Al)
    clear a;
    clear b;
    [a b]=xlsread('c:\rAlignmentBasePairComparison 1U6BCHAINB 1Y0QCHAINA',Al(k).SheetName);
   for i=1:length(b(:,2))
      b(i,2)=lower(strtrim(b(i,2)));
      b(i,5)=lower(strtrim(b(i,5)));
   end
   ct=0;
   for i=1:length(b(:,2))
      if isequal(b{i,2},'tsh')
         b{i,2}='ths';
      elseif isequal(b{i,2},'csh')
         b{i,2}='chs';
      elseif isequal(b{i,2},'tsw')
         b{i,2}='tws';
      elseif isequal(b{i,2},'csw')
         b{i,2}='cws';
      elseif isequal(b{i,2},'thw')
         b{i,2}='twh';
      elseif isequal(b{i,2},'chw')
         b{i,2}='cwh';
      elseif isequal(b{i,2},'ntsh')
         b{i,2}='nths';
      elseif isequal(b{i,2},'ncsh')
         b{i,2}='nchs';
      elseif isequal(b{i,2},'ntsw')
         b{i,2}='ntws';
      elseif isequal(b{i,2},'ncsw')
         b{i,2}='ncws';
      elseif isequal(b{i,2},'nthw')
         b{i,2}='ntwh';
      elseif isequal(b{i,2},'nchw')
         b{i,2}='ncwh';
      elseif ~isempty(b{i,2}) && isequal(b{i,2}(1),'s')
         b{i,2}='';
      end
   end
   for i=1:length(b(:,5))
      if isequal(b{i,5},'tsh')
         b{i,5}='ths';
      elseif isequal(b{i,5},'csh')
         b{i,5}='chs';
      elseif isequal(b{i,5},'tsw')
         b{i,5}='tws';
      elseif isequal(b{i,5},'csw')
         b{i,5}='cws';
      elseif isequal(b{i,5},'thw')
         b{i,5}='twh';
      elseif isequal(b{i,5},'chw')
         b{i,5}='cwh';
      elseif isequal(b{i,5},'ntsh')
         b{i,5}='nths';
      elseif isequal(b{i,5},'ncsh')
         b{i,5}='nchs';
      elseif isequal(b{i,5},'ntsw')
         b{i,5}='ntws';
      elseif isequal(b{i,5},'ncsw')
         b{i,5}='ncws';
      elseif isequal(b{i,5},'nthw')
         b{i,5}='ntwh';
      elseif isequal(b{i,5},'nchw')
         b{i,5}='ncwh';
      elseif ~isempty(b{i,5}) && isequal(b{i,5}(1),'s')
         b{i,5}='';
      end
   end

   numbp1=0;
   numnbp1=0;
   numnothing1=0;
   for i=1:length(b(:,2))
      len=length(b{i,2});
      if len==3
         numbp1=numbp1+1;
      elseif len==4 && ~isequal(b{i,2},'----')
         numnbp1=numnbp1+1;
      else 
         numnothing1=numnothing1+1;
      end
   end
%    [numbp1 numnbp1 numnothing1]

   numbp2=0;
   numnbp2=0;
   numnothing2=0;
   for i=1:length(b(:,5))
      len=length(b{i,5});
      if len==3
         numbp2=numbp2+1;
      elseif len==4 && ~isequal(b{i,5},'----')
         numnbp2=numnbp2+1;
      else 
         numnothing2=numnothing2+1;
      end
   end
%   [numbp2 numnbp2 numnothing2]

%We're comparing columns two and five.  
% Different Scenarios:
% 1) same basepairs
% 2) Different basepairs
% 3) Basepair in A, correct near bp in B
% 4) Bp in A, different near bp in B
% 5) Bp in A, aligned with 2 nt's making no bp in B
% 6) Bp in A, neither is aligned
% 7) Bp in A, only one is aligned
% 8) Basepair in B, correct near bp in A
% 9) Bp in B, different near bp in B
% 10) Bp in B, aligned with 2 nt's making no bp in A 
% 11) Bp in B, neither is aligned
% 12) Bp in B, only one is aligned
% 13) Single nucleotide in A aligned with nothing
% 14) Single nucleotide in B aligned with nothing
% 15) Single nucleotides aligned
% 16) Nbp in A with Nbp in B
% 17) Nbp in A, not with Nbp in B
% 18) Nbp in B, not with Nbp in A
   samebp = 0;
   diffbp = 0;
   inAnotB = 0;
   AwithNear = 0;
   AwithWrongNear = 0;
   AwithNoNT = 0;
   AwithOneNT = 0;
   inBnotA = 0;
   BwithNear = 0;
   BwithWrongNear = 0;
   BwithNoNT = 0;
   BwithOneNT = 0;
   SingleAwithNot = 0;
   SingleBwithNot = 0;
   SinglesAligned = 0;
   twoNears = 0;
   NearInAnotWithNearinB = 0;
   NearInBnotWithNearinA = 0;
   for i=1:length(b(:,2))
      if length(b{i,2})==3    
         if length(b{i,5})==3
            if isequal(b{i,2},b{i,5})
               samebp=samebp+1;                           % Case 1
            else
               diffbp=diffbp+1;                           % Case 2
            end
         elseif length(b{i,5})==4
            if isequal(b{i,5},'----')
               inAnotB=inAnotB+1;                        % Case 5
            elseif isequal(strcat('n',b{i,2}),b{i,5}) 
               AwithNear=AwithNear+1;                    % Case 3
            else
               AwithWrongNear=AwithWrongNear+1;          % Case 4
            end
         elseif isequal(b{i,5},'')
            if isequal(b{i,4},'---') && isequal(b{i,6},'---')
               AwithNoNT = AwithNoNT + 1;                % Case 6
            elseif isequal(b{i,4},'---') || isequal(b{i,6},'---')
               AwithOneNT = AwithOneNT + 1;              % Case 7
            end
         end
      elseif length(b{i,5})==3
         if length(b{i,2})==4
            if isequal(b{i,2},'----')
               inBnotA=inBnotA+1;                        % Case 10
            elseif isequal(strcat('n',b{i,5}),b{i,2})
               BwithNear=BwithNear+1;                    % Case 8
            else
               BwithWrongNear=BwithWrongNear+1;          % Case 9
            end
         elseif isequal(b{i,2},'')
            if isequal(b{i,1},'---') && isequal(b{i,3},'---')
               BwithNoNT = BwithNoNT + 1;                % Case 11
            elseif isequal(b{i,1},'---') || isequal(b{i,3},'---')
               BwithOneNT = BwithOneNT + 1;              % Case 12
            end
         end
      elseif isequal(b{i,4},'---') && isequal(b{i,6},'---')
         SingleAwithNot = SingleAwithNot + 1;            % Case 13
      elseif isequal(b{i,1},'---') && isequal(b{i,3},'---')
         SingleBwithNot = SingleBwithNot + 1;            % Case 14
      elseif isequal(b{i,3},'---') && isequal(b{i,6},'---') && ~isequal(b{i,1},'---') && ~isequal(b{i,4},'---')
         SinglesAligned = SinglesAligned + 1;            % Case 15
      elseif ~isempty(b{i,2}) && isequal(b{i,2}(1),'n')
         if isempty(b{i,5})     
             NearInAnotWithNearinB = NearInAnotWithNearinB + 1;
         elseif isequal(b{i,5}(1),'n')
             twoNears = twoNears+1;
         else
             NearInAnotWithNearinB = NearInAnotWithNearinB + 1;
         end
      elseif isequal(b{i,5}(1),'n')
         NearInBnotWithNearinA = NearInBnotWithNearinA + 1;
      end
   end
   T{2,k+1} = Al(k).Name;
   T{3,k+1} = samebp;
   T{4,k+1} = diffbp;
   T{5,k+1} = inAnotB;
   T{6,k+1} = AwithNear;
   T{7,k+1} = AwithWrongNear;
   T{8,k+1} = AwithNoNT;
   T{9,k+1} = AwithOneNT;
   T{10,k+1} = inBnotA;
   T{11,k+1} = BwithNear;
   T{12,k+1} = BwithWrongNear;
   T{13,k+1} = BwithNoNT;
   T{14,k+1} = BwithOneNT;
   T{15,k+1} = SingleAwithNot;
   T{16,k+1} = SingleBwithNot;
   T{17,k+1} = SinglesAligned;
   T{18,k+1} = twoNears;
   T{19,k+1} = NearInAnotWithNearinB;
   T{20,k+1} = NearInBnotWithNearinA;
end

% samebp
% diffbp
% inAnotB
% AwithNear
% AwithWrongNear
% AwithNoNT
% AwithOneNT
% inBnotA
% BwithNear
% BwithWrongNear
% BwithNoNT
% BwithOneNT
% SingleAwithNot
% SingleBwithNot
% SinglesAligned
File1=G(1);
File2=G(2);
T{2,1}=['Data for Alignment of ' File1.Filename ' and ' File2.Filename];
T{3,1}='Matching Basepairs Aligned';
T{4,1}='Different Basepairs Aligned';
T{5,1}=['Basepairs in ' File1.Filename ' with 2 NTs in ' File2.Filename ' not interacting'];
T{6,1}=['Basepairs in ' File1.Filename ' with matching near-bp in ' File2.Filename];
T{7,1}=['Basepairs in ' File1.Filename ' with non-matching near-bp in ' File2.Filename];
T{8,1}=['Basepairs in ' File1.Filename ' aligned with nothing in ' File2.Filename];
T{9,1}=['Basepairs in ' File1.Filename ' with one NT in ' File2.Filename];
T{10,1}=['Basepairs in ' File2.Filename ' with 2 NTs in ' File1.Filename ' not interacting'];
T{11,1}=['Basepairs in ' File2.Filename ' with matching near-bp in ' File1.Filename];
T{12,1}=['Basepairs in ' File2.Filename ' with non-matching near-bp in ' File1.Filename];
T{13,1}=['Basepairs in ' File2.Filename ' aligned with nothing in ' File1.Filename];
T{14,1}=['Basepairs in ' File2.Filename ' with one NT in ' File1.Filename];
T{15,1}=['Single NT in ' File1.Filename ' aligned with nothing in ' File2.Filename];
T{16,1}=['Single NT in ' File2.Filename ' aligned with nothing in ' File1.Filename];
T{17,1}=['Single NT in ' File1.Filename ' aligned with single NT in ' File2.Filename];
T{18,1}=['Near BP in ' File1.Filename ' near BP in ' File2.Filename];
T{19,1}=['Near BP in ' File1.Filename ' not aligned with near BP in ' File2.Filename];
T{20,1}=['Near BP in ' File2.Filename ' not aligned with near BP in ' File1.Filename];

xlswrite('Intron_alignment_data_2.xls',T);

% bothblank+match+nearmatch+notmatch
% samebp+diffbp+inAnotB+AwithNear+AwithWrongNear+AwithNoNT+AwithOneNT+inBnotA...
% +BwithNear+BwithWrongNear+BwithNoNT+BwithOneNT+SingleAwithNot+SingleBwithNot+SinglesAligned;
