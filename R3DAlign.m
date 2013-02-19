% File1 and File2 and can pdb id's such as '1s72' and '2j01' or they can be
% the variables containing the molecule information returned from
% zAddNTData.
% Chain1 is a cell array containing the chains (in character format) of the
%   fragments to be aligned.  Examples: {'1'},{'A'},{'3','1'}
% NTList1 is a cell the same size as Chain1 an
%   contains the nucleotide numbers of the fragments.  If all nucleotides
%   from a chain are to be aligned, {'all'} is an acceptable input.  Other
%   examples: {'all','all'} (if 2 chains are specified); and, {'2:5,7,8'}
% discCut is the discrepancy cutoff; use a cell array for iterative alignment
% numNeigh is the number of neighborhoods for each nucleotide
% cliqueMethod is either 'greedy' or 'full'; use a cell array for iteration
% Query.Type is either 'web' or 'local'.  Query.Name is the 13 character id
%   set by WebR3DAlign
% seed1 and seed2 are optional seed alignments. They can be ...

function [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,seed1,seed2)
%
tStart=tic;
if nargin < 11
    Query.Type = 'local';
end
try
if ~isfield(Query,'currIter')
   if iscell(discCut)
      if ~exist('seed1') || isempty(seed1)
         [AlignedNTs1,AlignedNTs2,ErrorMsg] = IterativeAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query);
      else
         [AlignedNTs1,AlignedNTs2,ErrorMsg] = IterativeAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,seed1,seed2);
      end
      return;
   else
      Query.currIter = 1; %Iteration Number
      discCut={discCut};
      numNeigh={numNeigh};
      bandwidth={bandwidth};
      cliqueMethod={cliqueMethod};
   end
else
   if Query.currIter > 1
       if ~exist('seed1') || isempty(seed1)
          exception = MException('MATLAB:nomem','Previous iteration had memory error');
          Indices1=[];
          Indices2=[];
          NTList{1}=NTList1;
          NTList{2}=NTList2;
          if ischar(File1),
             OutFilename = File1;
          else
             OutFilename = upper(File1.Filename);
          end
          for i=1:length(NTList1)
             OutFilename = [OutFilename '(' Chain1{i} ')' NTList1{i}]; %#ok<AGROW> 
          end
          if ischar(File2),
             OutFilename = [OutFilename '--' File2];
          else
             OutFilename = [OutFilename '--' File2.Filename];
          end
          for i=1:length(NTList2)
             OutFilename = [OutFilename '(' Chain2{i} ')' NTList2{i}]; %#ok<AGROW>
          end
          OutFilename=strrep(OutFilename, ':', '-');
          ShortOutFilename=OutFilename;
          for i=1:Query.currIter
             OutFilename = [OutFilename '_d' num2str(discCut{i}) '_p' num2str(numNeigh{i}) '_B' num2str(bandwidth{i})]; %#ok<AGROW>
          end
          OutFilename=strrep(OutFilename, '.', '');
          throw(exception);
       end
   end
end
addpath(genpath([pwd filesep 'FR3D']));
addpath([pwd filesep 'R3DAlign']);
if ~(exist([pwd filesep 'PDBFiles']) == 7),        % if directory doesn't yet exist
   mkdir([pwd filesep 'PDBFiles']);
end
addpath([pwd filesep 'PDBFiles']);
if ~(exist([pwd filesep 'PrecomputedData']) == 7),        % if directory doesn't yet exist
   mkdir([pwd filesep 'PrecomputedData']);
end
addpath([pwd filesep 'PrecomputedData']);
if ~(exist([pwd filesep 'Neighborhoods']) == 7),        % if directory doesn't yet exist
   mkdir([pwd filesep 'Neighborhoods']);
end
if ~(exist([pwd filesep 'Sequence Alignments']) == 7),        % if directory doesn't yet exist
   mkdir([pwd filesep 'Sequence Alignments']);
end
if ~(exist([pwd filesep 'R3D Align Output']) == 7),        % if directory doesn't yet exist
   mkdir([pwd filesep 'R3D Align Output']);
end
if ~(exist(fullfile(pwd, 'R3D Align Output', 'Bar Diagrams')) == 7),        % if directory doesn't yet exist
   mkdir(fullfile(pwd, 'R3D Align Output', 'Bar Diagrams'));
end
if ~(exist(fullfile(pwd, 'R3D Align Output', 'Spreadsheets')) == 7),        % if directory doesn't yet exist
   mkdir(fullfile(pwd, 'R3D Align Output', 'Spreadsheets'));
end
if ~(exist(fullfile(pwd, 'R3D Align Output', 'Text Alignments')) == 7),        % if directory doesn't yet exist
   mkdir(fullfile(pwd, 'R3D Align Output', 'Text Alignments'));
end
if ~(exist(fullfile(pwd, 'R3D Align Output', 'Final Mat Files')) == 7),        % if directory doesn't yet exist
   mkdir(fullfile(pwd, 'R3D Align Output', 'Final Mat Files'));
end
if ~(exist(fullfile(pwd, 'R3D Align Output', 'Summary Spreadsheets')) == 7),        % if directory doesn't yet exist
   mkdir(fullfile(pwd, 'R3D Align Output', 'Summary Spreadsheets'));
end
if nargin < 5
   discCut{Query.currIter} = input('Enter the discrepancy cutoff value: ');
end
if nargin < 6
   numNeigh{Query.currIter} = input('Enter the number of neighborhoods to retain for each nucleotide: ');
end
if nargin < 7
   bandwidth{Query.currIter} = input('Enter the band width for the seed alignment: ');
end
if nargin < 8
   cliqueMethod{Query.currIter} = input('Enter final clique method (Full or Greedy): ','s');
end

ErrorMsg='';
AlignedNTs1=cell(2);
AlignedNTs2=cell(2);

NTList{1}=NTList1;
NTList{2}=NTList2;

try
   if ischar(File1),
     fprintf('Loading PDB Info...\n');
     if isequal(File1,'uploaded')
        Filename1 = [Query.Name '_1'];
        File1 = zAddNTData(Filename1,0);
        File1.Filename=Query.UploadName1;
     else
        Filename1 = upper(File1);
        File1 = zAddNTData(Filename1,0);
     end
   else
      Filename1 = upper(File1.Filename);
   end
catch ME1
   ME1.identifier
   ME1.message
   ErrorMsg='ERROR: Unable to load PDB file for Molecule 1.';
   if isequal(Query.Type,'web')
      save([pwd filesep Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
   end
   return;
end

if File1.NumNT == -1
   ErrorMsg='ERROR: Unable to load PDB file for Molecule 1.';
   if isequal(Query.Type,'web')
      save([pwd filesep Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
   end
   return;
end

try
   if ischar(File2),
      if isequal(File2,'uploaded')
         Filename2 = [Query.Name '_2'];
         File2 = zAddNTData(Filename2,0);
         File2.Filename=Query.UploadName2;
      else
         Filename2 = upper(File2);
         File2 = zAddNTData(Filename2,0);
      end
   else
      Filename2 = upper(File2.Filename);
   end
catch ME2
   ME2.identifier
   ME2.message
   ErrorMsg='ERROR: Unable to load PDB file for Molecule 2.';
   if isequal(Query.Type,'web')
      save([pwd filesep Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
   end
   return;
end

if File1.NumNT == -1
   ErrorMsg='ERROR: Unable to load PDB file for Molecule 2.';
   if isequal(Query.Type,'web')
     save([pwd filesep Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
   end
   return;
end

Indices1=[];
Indices2=[];

OutFilename = Filename1;
NeighAFilename = ['Neighborhoods' filesep Filename1];
for i=1:length(NTList1)
   OutFilename = [OutFilename '(' Chain1{i} ')' NTList1{i}]; %#ok<AGROW>
   NeighAFilename = [NeighAFilename '(' Chain1{i} ')' NTList1{i}]; %#ok<AGROW>
   if isequal(Chain1{i},'all')
      Indices1=1:File1.NumNT;
   elseif isequal(NTList1{i},'all')
      chain=Chain1{i};
      c=cat(2,File1.NT.Chain);
      tmpIndices = find(lower(c)==lower(chain));
      Indices1 = [Indices1 tmpIndices]; %#ok<AGROW>
   else
      tmpIndices = zIndexLookup(File1,NTList1{i},Chain1{i});
      Indices1 = [Indices1 tmpIndices]; %#ok<AGROW>
   end
end
NeighAFilename = [NeighAFilename '_' num2str(numNeigh{Query.currIter}) '.mat'];
NeighAFilename = strrep(NeighAFilename, ':', '-');
NeighAFilename = [pwd filesep NeighAFilename];
OutFilename = [OutFilename '--' Filename2];
NeighBFilename = ['Neighborhoods' filesep Filename2];
for i=1:length(NTList2)
   OutFilename = [OutFilename '(' Chain2{i} ')' NTList2{i}]; %#ok<AGROW>
   NeighBFilename = [NeighBFilename '(' Chain2{i} ')' NTList2{i}]; %#ok<AGROW>
   if isequal(Chain2{i},'all')
      Indices2=1:File2.NumNT;
   elseif isequal(NTList2{i},'all')
      chain=Chain2{i};
      c=cat(2,File2.NT.Chain);
      tmpIndices = find(lower(c)==lower(chain));
      Indices2 = [Indices2 tmpIndices]; %#ok<AGROW>
   else
      tmpIndices = zIndexLookup(File2,NTList2{i},Chain2{i});
      Indices2 = [Indices2 tmpIndices]; %#ok<AGROW>
   end
end
OutFilename = strrep(OutFilename, ':', '-');
ShortOutFilename=OutFilename;
for i=1:Query.currIter
   OutFilename = [OutFilename '_d' num2str(discCut{i}) '_p' num2str(numNeigh{i}) '_B' num2str(bandwidth{i})]; %#ok<AGROW>
end
OutFilename = strrep(OutFilename, '.', '');

NeighBFilename = [NeighBFilename '_' num2str(numNeigh{Query.currIter}) '.mat'];
NeighBFilename = strrep(NeighBFilename, ':', '-');
NeighBFilename = [pwd filesep NeighBFilename];
if exist(fullfile(pwd, 'R3D Align Output','Final Mat Files', [OutFilename '.mat']))==2 %#ok<EXIST>
   disp('loading Final Alignment')
   load(fullfile(pwd, 'R3D Align Output', 'Final Mat Files', [OutFilename '.mat']));
else
   disp('not loading Final Alignment')
   if numNeigh{Query.currIter} > nchoosek(length(Indices1)-1,3)
      ErrorMsg='ERROR: Neighborhood parameter too high based on size of structure.';
      if isequal(Query.Type,'web')
         save([pwd Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
      end
      AlignedNTs1=[];
      AlignedNTs2=[];
      return;
   elseif numNeigh{Query.currIter} > nchoosek(length(Indices2)-1,3)
      ErrorMsg='ERROR: Neighborhood parameter too high based on size of structure.';
      if isequal(Query.Type,'web')
         save([pwd filesep Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
      end
      AlignedNTs1=[];
      AlignedNTs2=[];
      return;
   end

   disp('get neighborhoods')

   if exist(NeighAFilename)==2
      disp('loading NeighAFilename')
      load(NeighAFilename);
      AQuads=Quads; %#ok<NODEF>
      clear Quads;
      c = cat(1,File1.NT(Indices1).Center);           % will need matrix A later
      try
         File1.Distance = full(zMutualDistance(c,Inf));
      catch Err
         disp('Error in zMutualDistance')
         Err.identifier
         Err.message
         return;
      end
      A = triu(File1.Distance);
   else
      maxdist=15;
      getMoreNeigh = true;
      while getMoreNeigh == true;
         disp('not loading NeighAFilename')
         c = cat(1,File1.NT(Indices1).Center);           % nucleotide centers
         File1.Distance = full(zMutualDistance(c,Inf));  % Distance matrix for A
         A = triu(File1.Distance);                       % Distance matrix is symmetrical;
         [rA cA] = find(A<maxdist & A>0);
         [S1 IX] = sort(rA);
         rA = S1;
         cA = cA(IX);
         length(rA);
         ATrips = DoublesToTriples([rA cA]);
         AQuads = TriplesToQuads(ATrips);
         clear rA; clear cA; clear ATrips;
         if isempty(AQuads)
            maxdist=maxdist+10;
            getMoreNeigh = true;
            clear AQuads;
         else
            getMoreNeigh = false;
         end
         AQuads = rGetQuads(AQuads,A,numNeigh{Query.currIter});
         Quads=AQuads; 
         save(NeighAFilename, 'Quads')
         clear Quads;
      end
   end

   if exist(NeighBFilename)==2
      disp('loading NeighBFilename')
      load(NeighBFilename);
      BQuads=Quads;
      clear Quads;
      d = cat(1,File2.NT(Indices2).Center);           % will need matrix B later
      File2.Distance = full(zMutualDistance(d,Inf));
      B = triu(File2.Distance);
   else
      maxdist=15;
      getMoreNeigh = true;
      while getMoreNeigh == true;
         disp('not loading NeighBFilename')
         d = cat(1,File2.NT(Indices2).Center);           % nucleotide centers
         File2.Distance = full(zMutualDistance(d,Inf));  % Distance matrix for B
         B = triu(File2.Distance);                       % don't want duplicate entries returned in next lines
         [rB cB] = find(B<maxdist & B>0);
         [S2 IX] = sort(rB);
         rB = S2;
         cB = cB(IX);
         BTrips = DoublesToTriples([rB cB]);
         BQuads = TriplesToQuads(BTrips);
         clear rB; clear cB; clear BTrips;
         if isempty(BQuads)
            maxdist=maxdist+10;
            getMoreNeigh = true;
            clear BQuads;
         else
            getMoreNeigh = false;
         end
         BQuads = rGetQuads(BQuads,B,numNeigh{Query.currIter});
         Quads=BQuads; %#ok<NASGU>
         save(NeighBFilename, 'Quads')
         clear Quads;
      end
   end

   if isequal(Query.Type,'web') && ~isequal(Query.SeedName,'')
      FASTA = zReadFASTA(['seed_upload_file' '.fasta']);
      delete(['./SeedAlignments/' Query.Name '.fasta']);
      File1Sequence=cat(2,File1.NT(Indices1).Base);
      File2Sequence=cat(2,File2.NT(Indices2).Base);

      if abs(length(FASTA(1).Sequence)-length(File1Sequence))>0 || abs(length(FASTA(2).Sequence)-length(File2Sequence))>0
         ErrorMsg='ERROR: Seed alignment does not correspond with pdb files.';
         save([pwd Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
         return;
      else
         align1=[];
         align2=[];
         ct1=0;
         ct2=0;
         for i=1:length(FASTA(1).Aligned)
            add=1;
            if FASTA(1).Aligned(i)~='-' && FASTA(1).Aligned(i)~='.'
               ct1=ct1+1;
            else
               add=0;
            end
            if FASTA(2).Aligned(i)~='-' && FASTA(2).Aligned(i)~='.'
               ct2=ct2+1;
            else
               add=0;
            end
            if add==1
               align1(end+1)=ct1; %#ok<AGROW>
               align2(end+1)=ct2; %#ok<AGROW>
            end
         end
      end
   elseif ~exist('seed1') || isempty(seed1) %#ok<EXIST>
      if exist(fullfile(pwd, 'Sequence Alignments', [ShortOutFilename '.mat']))==2 %#ok<EXIST>
         disp('loading Sequence Alignment')
         load(fullfile(pwd, 'Sequence Alignments', [ShortOutFilename '.mat']));
      else
         disp('not loading Sequence Alignment')
         [align1 align2 charAlign1 charAlign2] = rGapNW(File1,Indices1,File2,Indices2,.999,7,0.25); %#ok<NASGU,ASGLU>
         save(fullfile(pwd, 'Sequence Alignments', [ShortOutFilename '.mat']), 'align1', 'align2', 'charAlign1', 'charAlign2')
      end
   else
      for i=1:length(seed1)
         if iscell(seed1)
            tmp=zIndexLookup(File1,seed1{i,1},seed1{i,2},0);
            align1(i)=find(Indices1==tmp); %#ok<AGROW>
            tmp=zIndexLookup(File2,seed2{i,1},seed2{i,2},0);
            align2(i)=find(Indices2==tmp); %#ok<AGROW>
         elseif isvector(seed1)
            align1=seed1;
            align2=seed2;
         end
      end
   end

%Files for the sequence alignment are created and moved to the
%'Sequence Alignments' folder.

   if ~isequal(Query.Type,'web') && exist(fullfile(pwd, 'Sequence Alignments',ShortOutFilename)) ~= 7
      I1=Indices1(align1);
      I2=Indices2(align2);
      rAlignmentSpreadsheet(File1,Indices1,File2,Indices2,I1,I2,ShortOutFilename,ErrorMsg);
      movefile([pwd filesep ShortOutFilename '*'], fullfile(pwd, 'Sequence Alignments', ShortOutFilename));

      clf
      [AAA,BBB] = rBarDiagram(File1,Indices1,File2,Indices2,I1,I2,ShortOutFilename,'R3D Align');
      View = [1 1 1 1 0 0 0];
      m1 = zBarDiagramInteractions(File1,Indices1,AAA,View,'above');
      m2 = zBarDiagramInteractions(File2,Indices2,BBB,View,'below');
      axis([0 20 m2(3) m1(4)])
      saveas(gcf,[ShortOutFilename '_int'],'pdf')

      movefile([pwd filesep ShortOutFilename '*'], fullfile(pwd, 'Sequence Alignments', ShortOutFilename));
      rWriteAlignmentFasta(File1,Indices1,File2,Indices2,I1,I2,NTList,ShortOutFilename);
      movefile([pwd filesep ShortOutFilename '*'], fullfile(pwd, 'Sequence Alignments', ShortOutFilename));
      rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,I1,I2,NTList,ShortOutFilename);
      movefile([pwd filesep ShortOutFilename '*'], fullfile(pwd, 'Sequence Alignments', ShortOutFilename));
      clear I1 I2;
   end

% Augment alignment with near matches
% For each nt in A aligned with a gap, the nearest aligned nt in B is found
% to use as the center of the band for that nt
  [align1,align2] = rAugmentAlignment(Indices1,align1,align2);

   SB=[];
   VMI=[];         %This will be a matrix containing the vertices that are created
                  %eg. if VMI(1) =[2 3], there was a vertex created for the
                  %2nd A neighborhood and 3rd B neighborhood
   nA=length(AQuads(:,1));
   nB=length(BQuads(:,1));

   List=cell(length(A),length(B));     %List{i,j} is a vector indicating the vertices that align
                                    %nucleotide i of A with nucleotide j of B

% Cutoff values for screening criterion
   cut2=4*4*2*discCut{Query.currIter}*discCut{Query.currIter};
   cut3=4*4*3*discCut{Query.currIter}*discCut{Query.currIter};
   cut4=4*4*4*discCut{Query.currIter}*discCut{Query.currIter};
   SA=round(bandwidth{Query.currIter}/2);
   fprintf('Making Local Alignment Graph...\n');
% Determine locations of vertices
   for i = 1:nA
      SB(1)=align2(AQuads(i,1));
      SB(2)=align2(AQuads(i,2));
      SB(3)=align2(AQuads(i,3));
      SB(4)=align2(AQuads(i,4));
      for j = 1:nB
         if abs(SB(1)-BQuads(j,1))<SA && abs(SB(2)-BQuads(j,2))<SA  && abs(SB(3)-BQuads(j,3))<SA && abs(SB(4)-BQuads(j,4))<SA
            screened=false;
            for k1 = 1:3
                if screened==false
                   for k2 = k1+1:4
                      AD=A(AQuads(i,[k1 k2]),AQuads(i,[k1 k2]));
                      BD=B(BQuads(j,[k1 k2]),BQuads(j,[k1 k2]));
                      DD=AD-BD;
                      DD=DD.*DD;
                      S=sum(sum(DD));
                      if S > cut2
                         screened = true;
                         break;
                      end
                   end
                end
            end

            if screened == false
               for k1 = 1:2
                  if screened == false
                     for k2 = k1+1:3
                        if screened == false
                           for k3 = k2+1:4
                              AD=A(AQuads(i,[k1 k2 k3]),AQuads(i,[k1 k2 k3]));
                              BD=B(BQuads(j,[k1 k2 k3]),BQuads(j,[k1 k2 k3]));
                              DD=AD-BD;
                              DD=DD.*DD;
                              S=sum(sum(DD));
                              if S > cut3
                                 screened=true;
                                 break;
                              end
                           end
                        end
                     end
                  end
               end
            end
            if screened == false
               AD=A(AQuads(i,:),AQuads(i,:));
               BD=B(BQuads(j,:),BQuads(j,:));
               DD=AD-BD;
               DD=DD.*DD;
               S=sum(sum(DD));
               if S > cut4
                  screened=true;
               end
            end
            if screened == false
               disc = xDiscrepancy(File1,Indices1(AQuads(i,:)),File2,Indices2(BQuads(j,:)));
               if disc < discCut{Query.currIter}
                  VMI=vertcat(VMI,[i j]); %#ok<AGROW>
                  List{AQuads(i,1),BQuads(j,1)}=[List{AQuads(i,1),BQuads(j,1)} length(VMI(:,1))];
                  List{AQuads(i,2),BQuads(j,2)}=[List{AQuads(i,2),BQuads(j,2)} length(VMI(:,1))];
                  List{AQuads(i,3),BQuads(j,3)}=[List{AQuads(i,3),BQuads(j,3)} length(VMI(:,1))];
                  List{AQuads(i,4),BQuads(j,4)}=[List{AQuads(i,4),BQuads(j,4)} length(VMI(:,1))];
               end
            end
         end
      end
   end

   if isempty(VMI)
      fprintf('\nNo nucleotides were aligned.\n');
      AlignedNTs1=[];
      AlignedNTs2=[];
      AlignedIndices1=[];
      AlignedIndices2=[];
      ErrorMsg='None aligned';
      if isequal(Query.Type,'web')
         save([pwd Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg', 'Query');
         rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
         rWriteAlignmentFasta(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
      else
         date = regexprep(datestr(now),':', '-');
         AlignmentName=['R3D Align ' File1.Filename ' ' File2.Filename ' ' date(1:17)];
         rWriteAlignmentFasta(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,AlignmentName);
         rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,AlignmentName);
      end
      return;
   end
   
try
   [EM]=rMakeEdgeMatrix(VMI,List);
catch Err
  disp('Error in rMakeEdgeMatrix')
  Err.identifier
  Err.message
%   Err.stack
  if (strcmp(Err.identifier,'MATLAB:nomem'))
     ErrorMsg='Out of Memory';
  end
  AlignedIndices1=[];
  AlignedIndices2=[];
  AlignedNTs1{1,1}=[];
  AlignedNTs2{1,1}=[];
  save(fullfile(pwd, 'R3D Align Output', 'Final Mat Files', [OutFilename '.mat']),'Indices1','Indices2','AlignedIndices1','AlignedIndices2');
  WriteOutput(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename,ShortOutFilename,ErrorMsg,Query)
  return;
end  
   clear List;
   clear A;
   clear B;

   maximalClique = rGetMaximalClique(EM);
   CLB = length(maximalClique);

   refineMore=true;
   numLoop=0;
   fprintf('Preprocessing...\n');
   while refineMore==true
      numLoop=numLoop+1;
      D=[];
      C=colorVertices(EM);       %vector of length=length(VMI) with color class for each vertex

      for p = 1:length(VMI(:,1))
         if ~any(maximalClique==p)
            CX = C(EM(p,:)==0);               %Colors of neighbors
            if length(unique(CX))<CLB-1
               D=[D p]; %#ok<AGROW>
            end
         end
      end
      if length(D)/length(VMI(:,1)) < .2 || numLoop == 2
         refineMore=false;
      end
      VMI(D,:)=[]; %#ok<AGROW>
      EM(D,:)=[];
      EM(:,D)=[];
   end
   fprintf('Finding Clique...\n');

   if strcmpi(cliqueMethod{Query.currIter},'full')
      MC = rCliqueBB(EM);
   elseif strcmpi(cliqueMethod{Query.currIter},'greedy')
      MC = rGetMaximalClique(EM);
   else
      method = input('Use full clique finding method? Y/N :','s');
      while ~strcmpi(method,'y') && ~strcmpi(method,'n')
         fprintf('Invalid user input.\n');
         method = input('Use full clique finding method? Y/N :','s');
      end
      if strcmpi(method,'y')
         MC = rCliqueBB(EM,'c');
      else
         MC = rGetMaximalClique(EM);
      end
   end

   fprintf('Producing Alignment Output...\n');
% Convert clique into corresponding alignment
   Z=VMI(MC,:);
   tmpAlmnt = [AQuads(Z(:,1),1) BQuads(Z(:,2),1);
               AQuads(Z(:,1),2) BQuads(Z(:,2),2);
               AQuads(Z(:,1),3) BQuads(Z(:,2),3);
               AQuads(Z(:,1),4) BQuads(Z(:,2),4)];  %unsorted, with duplicates
   [S1 IX] = sort(tmpAlmnt(:,1));
   tmpAlmnt(:,1)=S1;
   tmpAlmnt(:,2) = tmpAlmnt(IX,2);    %sorted
   [b m] = unique(tmpAlmnt(:,1));
   Indnums1=b;
   Indnums2=tmpAlmnt(m,2);            %no duplicates

%The following are the indices of aligned nucleotides
   AlignedIndices1 = Indices1(Indnums1);
   AlignedIndices2 = Indices2(Indnums2);

%AlignedNTs contains the nucleotide number, base and chain of each aligned nucleotide
   clear AlignedNTs1;
   clear AlignedNTs2;
   AlignedNTs1=cell(length(AlignedIndices1),2);
   AlignedNTs2=cell(length(AlignedIndices1),2);
   for i=1:length(AlignedIndices1)
      AlignedNTs1{i,1}=File1.NT(AlignedIndices1(i)).Number;
      AlignedNTs1{i,2}=File1.NT(AlignedIndices1(i)).Chain;
      AlignedNTs1{i,3}=File1.NT(AlignedIndices1(i)).Base;
      AlignedNTs1{i,4}=AlignedIndices1(i);
      
      AlignedNTs2{i,1}=File2.NT(AlignedIndices2(i)).Number;
      AlignedNTs2{i,2}=File2.NT(AlignedIndices2(i)).Chain;
      AlignedNTs2{i,3}=File2.NT(AlignedIndices2(i)).Base;
      AlignedNTs2{i,4}=AlignedIndices2(i);
   end

   save(fullfile(pwd, 'R3D Align Output', 'Final Mat Files', [OutFilename '.mat']),'Indices1','Indices2','AlignedIndices1','AlignedIndices2','AlignedNTs1','AlignedNTs2');
end
Query.Time=toc(tStart);
WriteOutput(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename,ShortOutFilename,ErrorMsg,Query)

catch ME
  ME.identifier
  ME.message
  AlignedIndices1=[];
  AlignedIndices2=[];
%   ME.stack
  if (strcmp(ME.identifier,'MATLAB:nomem'))
     Query.ErrorMsg = 'Out of Memory';
     ErrorMsg='Out of Memory';
     Query.Time=toc(tStart);
     AlignedNTs1{1,1}=[];
     AlignedNTs2{1,1}=[];
     save(fullfile(pwd, 'R3D Align Output', 'Final Mat Files', [OutFilename '.mat']),'Indices1','Indices2','AlignedIndices1','AlignedIndices2');
     WriteOutput(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename,ShortOutFilename,ErrorMsg,Query)
  else
     for k=1:length(ME.stack)
       ME.stack(k)
     end
  end
  ErrorMsg = ME.message;
  return;
end
end

function WriteOutput(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename,ShortOutFilename,ErrorMsg,Query)
if isequal(Query.Type,'web') 
   if ~strcmp(ErrorMsg,'Out of Memory')
      if length(AlignedIndices1)>4
         clf
         [AAA,BBB] = rBarDiagram(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,Query.Name,'R3D Align');
         View = [1 1 1 1 0 0 0];
         m1 = zBarDiagramInteractions(File1,Indices1,AAA,View,'above');
         m2 = zBarDiagramInteractions(File2,Indices2,BBB,View,'below');
         axis([0 20 m2(3) m1(4)])
         saveas(gcf,[Query.Name '_int'],'pdf')
         saveas(gcf,[Query.Name '_int'],'png')      
      end
      rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
      rWriteAlignmentFasta(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
      VP.Write=1;
      rSuperimposeNucleotides(File1,[AlignedIndices1 setdiff(Indices1,AlignedIndices1)],File2,[AlignedIndices2 setdiff(Indices2,AlignedIndices2)],VP,length(AlignedIndices1),Query.Name);
      FL=rAlignmentSpreadsheet(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,Query.Name,ErrorMsg);
   %    rWriteSummaryStatistics(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,Query.Name,FL);
      save([pwd filesep Query.Name '.mat'], 'AlignedIndices1', 'AlignedIndices2', 'ErrorMsg', 'Query');
   end
elseif exist(fullfile(pwd, 'R3D Align Output', OutFilename)) ~= 7 %if folder does not exist
   if length(AlignedIndices1)>4
      clf
      [AAA,BBB] = rBarDiagram(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,OutFilename,'R3D Align');
      View = [1 1 1 1 0 0 0];
      m1 = zBarDiagramInteractions(File1,Indices1,AAA,View,'above');
      m2 = zBarDiagramInteractions(File2,Indices2,BBB,View,'below');
      axis([0 20 m2(3) m1(4)])
      saveas(gcf,[OutFilename '_int'],'pdf')

      movefile([pwd filesep OutFilename '*'], fullfile(pwd, 'R3D Align Output', OutFilename));
%       copyfile(fullfile(pwd,'R3D Align Output',OutFilename,[OutFilename '.pdf']), fullfile(pwd, 'R3D Align Output', 'Bar Diagrams'));
   end
   rWriteAlignmentFasta(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename);
   movefile([pwd filesep OutFilename '*'], fullfile(pwd, 'R3D Align Output', OutFilename));
   FL=rAlignmentSpreadsheet(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,OutFilename,ErrorMsg);
   movefile([pwd filesep OutFilename '*'], fullfile(pwd, 'R3D Align Output', OutFilename));
   rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename);
   movefile([pwd filesep OutFilename '*'], fullfile(pwd, 'R3D Align Output', OutFilename));
%    copyfile(fullfile(pwd,'R3D Align Output',OutFilename,[OutFilename '.txt']), fullfile(pwd, 'R3D Align Output', 'Text Alignments'));   
   rAnalyzeAlignmentNew(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,OutFilename,ShortOutFilename,ErrorMsg,Query); 
%    rWriteSummaryStatistics(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,OutFilename,FL);
%    if exist(fullfile(pwd, 'R3D Align Output', OutFilename, [OutFilename '.xlsx'])) == 2
%       copyfile(fullfile(pwd,'R3D Align Output',OutFilename,[OutFilename '.xlsx']), fullfile(pwd, 'R3D Align Output', 'Spreadsheets'));
%    elseif exist(fullfile(pwd, 'R3D Align Output', OutFilename, [OutFilename '.xls'])) == 2
%       copyfile(fullfile(pwd,'R3D Align Output',OutFilename,[OutFilename '.xls']), fullfile(pwd, 'R3D Align Output', 'Spreadsheets'));
%    else
%       disp('File not found for rAnalyzeAlignment');
%    end
end                  
end

function [AlignedNTs1,AlignedNTs2,ErrorMsg] =IterativeAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,seed1,seed2)
   L=length(discCut);
   Query.currIter = 1; %current Iteration
   if ~exist('seed1') || isempty(seed1)
      [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query);
   else
      [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,seed1,seed2);
   end
   ct=2;
   while ct <= L
      Query.currIter = ct;
      if ~isempty(AlignedNTs1{1,1})
         [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,AlignedNTs1,AlignedNTs2);
      else
         [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query);
      end
      ct = ct + 1;
   end
end
