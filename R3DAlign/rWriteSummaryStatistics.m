% This file has been adopted from EvalConsensus.m and
% zCompare16SAlignments;
% It could be written more generally and efficiently so that it isn't based
% on the spreadsheet computed by rAlignmentSpreadsheet.
% Evaluate the number of near basepairs aligned with true basepairs, near with near, 
% true with true, etc.  Evaluate # of conserved nested cWW, non-nested,
% etc.  Evaluate # of nt's aligned, # of base matches.
function [] = rWriteSummaryStatistics(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,QueryName,FL)

   b=FL;
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


   samebp = 0;
   diffbp = 0;
   for i=1:length(b(:,2))
      if length(b{i,2})==3 && length(b{i,5})==3
            if isequal(b{i,2},b{i,5})
               samebp=samebp+1;                           % Case 1
            else
               diffbp=diffbp+1;                           % Case 2
            end
      end
   end
   
   numbp=samebp+diffbp;
   
   if length(AlignedIndices1) > 4
      Discrep = rFindAlignmentDiscrepancies(File1,AlignedIndices1,File2,AlignedIndices2,'nearest4');
   elseif length(AlignedIndices1) == 4
      Discrep = xDiscrepancy(File1,AlignedIndices1,File2,AlignedIndices2);
      Discrep = mean(Discrep);
   else
      Discrep = 'N/A';
   end
   
   if ~isempty(AlignedIndices1) 
      GlobDiscrep = xDiscrepancy(File1,AlignedIndices1,File2,AlignedIndices2);
   else
      GlobDiscrep = 'N/A';
   end
   
   Filename=[QueryName '_stats.csv'];
   fidOUT = fopen(Filename,'w+');
   fprintf(fidOUT,'Number of nucleotides in structure 1,%s\r\n', num2str(length(Indices1)));
   fprintf(fidOUT,'Number of nucleotides in structure 2,%s\r\n', num2str(length(Indices2)));
   fprintf(fidOUT,'Number of nucleotides aligned,%s\r\n', num2str(length(AlignedIndices1)));
   fprintf(fidOUT,'Percentage of structure 1 nucleotides aligned,%s\r\n', num2str(100*length(AlignedIndices1)/length(Indices1)));
   fprintf(fidOUT,'Percentage of structure 2 nucleotides aligned,%s\r\n', num2str(100*length(AlignedIndices2)/length(Indices2)));
   fprintf(fidOUT,'Number of basepairs in structure 1,%s\r\n', num2str(numbp1));
   fprintf(fidOUT,'Number of basepairs in structure 2,%s\r\n', num2str(numbp2));
   fprintf(fidOUT,'Number of basepairs aligned,%s\r\n', num2str(numbp));
   fprintf(fidOUT,'Percentage of structure 1 basepairs aligned,%s\r\n', num2str(numbp/numbp1*100));
   fprintf(fidOUT,'Percentage of structure 2 basepairs aligned,%s\r\n', num2str(numbp/numbp2*100));
   fprintf(fidOUT,'Mean local neighborhood discrepancy,%s\r\n', Discrep);
   fprintf(fidOUT,'Global discrepancy of all aligned nucleotides,%s\r\n', GlobDiscrep);
   fclose(fidOUT);