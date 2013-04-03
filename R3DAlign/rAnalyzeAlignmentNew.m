% This file has been adopted from EvalConsensus.m and
% zCompare16SAlignments;
% Evaluate the number of near basepairs aligned with true basepairs, near with near, 
% true with true, etc.  Evaluate # of conserved nested cWW, non-nested,
% etc.  Evaluate # of nt's aligned, # of base matches.
function [T] = rAnalyzeAlignmentNew(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,OutFilename,ShortOutFilename,ErrorMsg,Query)
if ~strcmp(ErrorMsg,'Out of Memory')  && ~strcmp(ErrorMsg,'None aligned') 
   if exist(fullfile(pwd, 'R3D Align Output', OutFilename, [OutFilename '.xlsx'])) == 2 %#ok<EXIST>
      SpreadsheetName = fullfile(pwd, 'R3D Align Output', OutFilename, [OutFilename '.xlsx']);
   elseif exist(fullfile(pwd, 'R3D Align Output', OutFilename, [OutFilename '.xls'])) == 2 %#ok<EXIST>
      SpreadsheetName = fullfile(pwd, 'R3D Align Output', OutFilename, [OutFilename '.xls']);
   else
      disp('File not found for rAnalyzeAlignment');
   end
   T=cell(29,2);
    clear a;
    clear b;
    [a b]=xlsread(SpreadsheetName); %#ok<ASGLU>
   for i=1:length(b(:,2))
      b(i,2)=lower(strtrim(b(i,2)));
      b(i,5)=lower(strtrim(b(i,5)));
   end
%    ct=0;
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

%------------------------------------------------------------------------
% --------------------------------------- Calculations for each alignment
  File(1)=File1;
  File(2)=File2;
  Tally = rTallyInteractions(File,AlignedIndices1,AlignedIndices2,0);
  
% --------------------------------------- Superimpose local neighborhoods
   if isempty(File(1).Distance),
      for f = 1:length(File),
         c = cat(1,File(f).NT(1:File(f).NumNT).Center); % nucleotide centers
         File(f).Distance = zMutualDistance(c,16); %#ok<AGROW> % compute distances < 16 Angstroms
      end
   end
  if length(AlignedIndices1) > 4
     Discrep = rFindAlignmentDiscrepancies(File1,AlignedIndices1,File2,AlignedIndices2,'nearest4');
  elseif length(AlignedIndices1) == 4
     Discrep = xDiscrepancy(File1,AlignedIndices1,File2,AlignedIndices2);
  end
%   Discrep = zHistogramDiscrepanciesInAlignment(File,AlignedIndices1,AlignedIndices2);


% -------------------------------- Calculate IDI of aligned base combinations 
%                                  with basepairs in 3D structure
%   IDI     = rAlignmentIsostericity(File,AlignedIndices1,AlignedIndices2,0);

  Identical = length(find(cat(1,File(1).NT(AlignedIndices1).Code)==cat(1,File(2).NT(AlignedIndices2).Code)));
  
   
  
   T{1,1}=['Data for Alignment of ' File1.Filename ' and ' File2.Filename];
   T{3,1} = 'Number aligned';
   T{4,1} = 'Number of exact base matches in alignment';
   T{5,1}='Basepairs From Same Geometeric Family Aligned';
   T{6,1}='Basepairs From Different Geometric Family Aligned';
   T{7,1} = 'Nested cWW aligned';
   T{8,1} = 'Nested non-cWW aligned';
   T{9,1} = 'Non-nested cWW aligned';
   T{10,1} = 'Non-nested non-cWW aligned';
   T{11,1} = 'Stacking aligned';
   T{12,1} = 'Base-phosphate aligned';
   T{13,1} = 'Mean local geometric discrepancy';
   T{14,1} = 'Median local geometric discrepancy';
   T{15,1} = 'Mean IDI between aligned base combinations and real pairs';
   T{16,1} = 'Median IDI between aligned base combinations and real pairs';
   T{17,1}=['Basepairs in ' File1.Filename ' with 2 NTs in ' File2.Filename ' not interacting'];
   T{18,1}=['Basepairs in ' File1.Filename ' with matching near-bp in ' File2.Filename];
   T{19,1}=['Basepairs in ' File1.Filename ' with non-matching near-bp in ' File2.Filename];
   T{20,1}=['Basepairs in ' File1.Filename ' aligned with nothing in ' File2.Filename];
   T{21,1}=['Basepairs in ' File1.Filename ' with one NT in ' File2.Filename];
   T{22,1}=['Basepairs in ' File2.Filename ' with 2 NTs in ' File1.Filename ' not interacting'];
   T{23,1}=['Basepairs in ' File2.Filename ' with matching near-bp in ' File1.Filename];
   T{24,1}=['Basepairs in ' File2.Filename ' with non-matching near-bp in ' File1.Filename];
   T{25,1}=['Basepairs in ' File2.Filename ' aligned with nothing in ' File1.Filename];
   T{26,1}=['Basepairs in ' File2.Filename ' with one NT in ' File1.Filename];
   T{27,1}=['Single NT in ' File1.Filename ' aligned with nothing in ' File2.Filename];
   T{28,1}=['Single NT in ' File2.Filename ' aligned with nothing in ' File1.Filename];
   T{29,1}=['Single NT in ' File1.Filename ' aligned with single NT in ' File2.Filename];
   % T{18,1}=['Near BP in ' File1.Filename ' near BP in ' File2.Filename];
   % T{19,1}=['Near BP in ' File1.Filename ' not aligned with near BP in ' File2.Filename];
   % T{20,1}=['Near BP in ' File2.Filename ' not aligned with near BP in ' File1.Filename];

   T{3,2} = length(AlignedIndices1);
   T{4,2} = Identical;
   T{5,2} = samebp/2;
   T{6,2} = diffbp/2;
   for v = 1:6,
     T{v+6,2} = Tally(v);
   end
   T{13,2} = mean(Discrep);
   T{14,2} = median(Discrep);
%    T{15,2} = mean(IDI);
%    T{16,2} = median(IDI);
   T{17,2} = inAnotB;
   T{18,2} = AwithNear;
   T{19,2} = AwithWrongNear;
   T{20,2} = AwithNoNT;
   T{21,2} = AwithOneNT;
   T{22,2} = inBnotA;
   T{23,2} = BwithNear;
   T{24,2} = BwithWrongNear;
   T{25,2} = BwithNoNT;
   T{26,2} = BwithOneNT;
   T{27,2} = SingleAwithNot;
   T{28,2} = SingleBwithNot;
   T{29,2} = SinglesAligned;
%    T{18,2} = twoNears;
%    T{19,2} = NearInAnotWithNearinB;
%    T{20,2} = NearInBnotWithNearinA;

if ~ismac
   xlswrite(SpreadsheetName,T,'Sheet1','I1');
   clear T;
end

end

if ~exist(fullfile(pwd,'R3D Align Output','Summary Spreadsheets',[ShortOutFilename '.xls']),'file') && ~exist(fullfile(pwd,'R3D Align Output','Summary Spreadsheets',[ShortOutFilename '.xls']),'file')
   T{1,1}='Molecule';
   T{2,1}=File1.Filename;
   T{3,1}=File2.Filename;
   T{1,2}='Number of Nucleotides';
   T{2,2}=length(Indices1);
   T{3,2}=length(Indices2);
   T{1,3}='Number of Basepairs';
   if ~strcmp(ErrorMsg,'Out of Memory') && ~strcmp(ErrorMsg,'None aligned')
      T{2,3}= samebp/2+diffbp/2+inAnotB+AwithNear+AwithWrongNear+AwithNoNT+AwithOneNT;
      T{3,3}= samebp/2+diffbp/2+inBnotA+BwithNear+BwithWrongNear+BwithNoNT+BwithOneNT;
   end
   T{5,2} = 'Number aligned';
   T{5,3} = 'Number of exact base matches in alignment';
   T{5,4}='Basepairs From Same Geometeric Family Aligned';
   T{5,5}='Basepairs From Different Geometric Family Aligned';
   T{5,6} = 'Nested cWW aligned';
   T{5,7} = 'Nested non-cWW aligned';
   T{5,8} = 'Non-nested cWW aligned';
   T{5,9} = 'Non-nested non-cWW aligned';
   T{5,10} = 'Stacking aligned';
   T{5,11} = 'Base-phosphate aligned';
   T{5,12} = 'Mean local geometric discrepancy';
   T{5,13} = 'Median local geometric discrepancy';
   T{5,14} = 'Mean IDI between aligned base combinations and real pairs';
   T{5,15} = 'Median IDI between aligned base combinations and real pairs';
   T{5,16}=['Basepairs in ' File1.Filename ' with 2 NTs in ' File2.Filename ' not interacting'];
   T{5,17}=['Basepairs in ' File1.Filename ' with matching near-bp in ' File2.Filename];
   T{5,18}=['Basepairs in ' File1.Filename ' with non-matching near-bp in ' File2.Filename];
   T{5,19}=['Basepairs in ' File1.Filename ' aligned with nothing in ' File2.Filename];
   T{5,20}=['Basepairs in ' File1.Filename ' with one NT in ' File2.Filename];
   T{5,21}=['Basepairs in ' File2.Filename ' with 2 NTs in ' File1.Filename ' not interacting'];
   T{5,22}=['Basepairs in ' File2.Filename ' with matching near-bp in ' File1.Filename];
   T{5,23}=['Basepairs in ' File2.Filename ' with non-matching near-bp in ' File1.Filename];
   T{5,24}=['Basepairs in ' File2.Filename ' aligned with nothing in ' File1.Filename];
   T{5,25}=['Basepairs in ' File2.Filename ' with one NT in ' File1.Filename];
   T{5,26}=['Single NT in ' File1.Filename ' aligned with nothing in ' File2.Filename];
   T{5,27}=['Single NT in ' File2.Filename ' aligned with nothing in ' File1.Filename];
   T{5,28}=['Single NT in ' File1.Filename ' aligned with single NT in ' File2.Filename];
   T{5,29}='Running Time';
  
   T{6,1} = OutFilename;
   if ~strcmp(ErrorMsg,'Out of Memory') && ~strcmp(ErrorMsg,'None aligned')
      T{6,2} = length(AlignedIndices1);
      T{6,3} = Identical;
      T{6,4} = samebp/2;
      T{6,5} = diffbp/2;
      for v = 1:6,
        T{6,v+5} = Tally(v);
      end
      T{6,12} = mean(Discrep);
      T{6,13} = median(Discrep);
%       T{6,14} = mean(IDI);
%       T{6,15} = median(IDI);
      T{6,16} = inAnotB;
      T{6,17} = AwithNear;
      T{6,18} = AwithWrongNear;
      T{6,19} = AwithNoNT;
      T{6,20} = AwithOneNT;
      T{6,21} = inBnotA;
      T{6,22} = BwithNear;
      T{6,23} = BwithWrongNear;
      T{6,24} = BwithNoNT;
      T{6,25} = BwithOneNT;
      T{6,26} = SingleAwithNot;
      T{6,27} = SingleBwithNot;
      T{6,28} = SinglesAligned;
      T{6,29} = Query.Time;
   else
      T{6,2} = ErrorMsg;
      for i=3:28
            T{6,i}=0;
      end
   end
      if ~ismac
         xlswrite(fullfile(pwd,'R3D Align Output','Summary Spreadsheets',[ShortOutFilename '.xls']),T,'Sheet1','A1');
      end
else
   if exist(fullfile(pwd,'R3D Align Output','Summary Spreadsheets',[ShortOutFilename '.xls']),'file')
      clear a b;
      [a b]=xlsread(fullfile(pwd,'R3D Align Output','Summary Spreadsheets',[ShortOutFilename '.xls'])); %#ok<ASGLU>
   else exist(fullfile(pwd,'R3D Align Output','Summary Spreadsheets',[ShortOutFilename '.xlsx']),'file')
      clear a b;
      [a b]=xlsread(fullfile(pwd,'R3D Align Output','Summary Spreadsheets',[ShortOutFilename '.xlsx'])); %#ok<ASGLU>
   end
   found=false;
   L=length(b(:,1));
   for i = 1:L
      if isequal(b{i,1},OutFilename)
         found = true;
         break;
      end
   end
   if ~found
      T{1,1} = OutFilename;
      if ~strcmp(ErrorMsg,'Out of Memory') && ~strcmp(ErrorMsg,'None aligned')
         T{1,2} = length(AlignedIndices1);
         T{1,3} = Identical;
         T{1,4} = samebp/2;
         T{1,5} = diffbp/2;
         for v = 1:6,
            T{1,v+5} = Tally(v);
         end
         T{1,12} = mean(Discrep);
         T{1,13} = median(Discrep);
%          T{1,14} = mean(IDI);
%          T{1,15} = median(IDI);
         T{1,16} = inAnotB;
         T{1,17} = AwithNear;
         T{1,18} = AwithWrongNear;
         T{1,19} = AwithNoNT;
         T{1,20} = AwithOneNT;
         T{1,21} = inBnotA;
         T{1,22} = BwithNear;
         T{1,23} = BwithWrongNear;
         T{1,24} = BwithNoNT;
         T{1,25} = BwithOneNT;
         T{1,26} = SingleAwithNot;
         T{1,27} = SingleBwithNot;
         T{1,28} = SinglesAligned;
         T{1,29} = Query.Time;
      else
         T{1,2}=ErrorMsg;
         for i=3:28
            T{1,i}=0;
         end
      end
      if ~ismac
         xlswrite(fullfile(pwd,'R3D Align Output','Summary Spreadsheets',[ShortOutFilename '.xls']),T,'Sheet1',['A' num2str(L+1)]);
      end
   end
end





