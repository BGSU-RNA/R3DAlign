% FileList is  a cell array of pdb ids
% e.g. Filelist = {'1s7s','2qbg','2j01'}

function [Good Bad] = CreateSeqAlignments(FileList)
% disp(num2str(length(FileList)))
% Good={};
% numGood=0;
% Bad={};
% numBad=0;
% for i1=1:length(FileList)
%    i1
%    Filename1 = upper(FileList{i1});  
%    File1 = zAddNTData(Filename1,0);
%    if ~isempty(File1.Backbone)
%        numGood = numGood + 1
%        Good{numGood} = Filename1; %#ok<AGROW>
%    else
%        numBad = numBad + 1
%        Bad{numBad} = Filename1; %#ok<AGROW>
%    end      
% end
% FileList = Good;
% disp(num2str(length(FileList)))

addpath(genpath([pwd filesep 'FR3D']));
addpath([pwd filesep 'R3DAlign']);
if ~(exist([pwd filesep 'PDBFiles']) == 7),        %#ok<*EXIST> % if directory doesn't yet exist
   mkdir([pwd filesep 'PDBFiles']);
end
addpath([pwd filesep 'PDBFiles']);
if ~(exist([pwd filesep 'PrecomputedData']) == 7),        % if directory doesn't yet exist
   mkdir([pwd filesep 'PrecomputedData']);
end
addpath([pwd filesep 'PrecomputedData']);
if ~(exist([pwd filesep 'Sequence Alignments']) == 7),        % if directory doesn't yet exist
   mkdir([pwd filesep 'Sequence Alignments']);
end

for i1=1:length(FileList)-1
   Filename1 = upper(FileList{i1});  
   File1 = zAddNTData(Filename1,0);  
   if ~isempty(File1.Backbone)
      Chains1 = unique(cat(2,File1.NT.Chain));
      for i2=i1+1:length(FileList)
         Filename2 = upper(FileList{i2});
         File2 = zAddNTData(Filename2,0); 
         if ~isempty(File2.Backbone)
            Chains2 = unique(cat(2,File2.NT.Chain)); 
            for j1=1:length(Chains1)
               c=cat(2,File1.NT.Chain);
               Indices1 = find(lower(c)==lower(Chains1(j1)));
               for j2=1:length(Chains2)
                  SeqFilename1 = [fullfile(pwd, 'Sequence Alignments', Filename1) '(' Chains1(j1) ')' 'all--' Filename2 '(' Chains2(j2) ')' 'all.mat'];
                  SeqFilename2 = [fullfile(pwd, 'Sequence Alignments', Filename2) '(' Chains2(j2) ')' 'all--' Filename1 '(' Chains1(j1) ')' 'all.mat'];
                  d=cat(2,File2.NT.Chain);
                  Indices2 = find(lower(d)==lower(Chains2(j2)));
                  if exist(SeqFilename1)==2 %#ok<EXIST>
                     disp([SeqFilename1 ' already exists.']);
                     if exist(SeqFilename2)~=2
                        disp(['Creating ' SeqFilename2]);
                        M = max(length(Indices1)/length(Indices2),length(Indices2)/length(Indices1));
                        if M < 1.5
                           if length(Indices1) > 300 || length(Indices2) > 300
                              [align1 align2 charAlign1 charAlign2] = rGapNW(File2,Indices2,File1,Indices1,.999,7,.25); %#ok<*NASGU,ASGLU>
                              save(SeqFilename2, 'align1', 'align2', 'charAlign1', 'charAlign2');
                           else
                              disp(['   Both sequences are less than 300 nucleotides. File not created']);
                           end
                        else
                           disp(['   Lengths of chains differ by a factor of ' num2str(M) '. File not created']);
                        end
                     else
                        disp([SeqFilename2 ' already exists.']); 
                     end
                  else
                     disp(['Creating ' SeqFilename1]);
                     M = max(length(Indices1)/length(Indices2),length(Indices2)/length(Indices1));
                     if M < 1.5
                        if length(Indices1) > 300 || length(Indices2) > 300
                           [align1 align2 charAlign1 charAlign2] = rGapNW(File2,Indices2,File1,Indices1,.999,7,.25); %#ok<*NASGU,ASGLU>
                           save(SeqFilename2, 'align1', 'align2', 'charAlign1', 'charAlign2');
                        else
                           disp(['   Both sequences are less than 300 nucleotides. File not created']);
                        end
                        if exist(SeqFilename2)==2
                           disp([SeqFilename2 ' already exists.']);
                        else
                           disp(['Creating ' SeqFilename2]);
                           tmpalign1=align1; tmpalign2=align2; tmpcharAlign1=charAlign1; tmpcharAlign2=charAlign2;
                           align1=tmpalign2; align2=tmpalign1; charAlign1=tmpcharAlign2; charAlign2=tmpcharAlign1;
                           save(SeqFilename2, 'align1', 'align2', 'charAlign1', 'charAlign2');
                        end
                     else
                        disp(['   Lengths of chains differ by a factor of ' num2str(M) '. File not created']);
                        disp(['Creating ' SeqFilename2]);
                        disp(['   Lengths of chains differ by a factor of ' num2str(M) '. File not created']);
                     end
                  end
               end
            end 
         end
      end
   end
end