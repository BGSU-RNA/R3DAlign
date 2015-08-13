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
% Query.Type is either 'web' or 'local'.  
% Query.Name is the 13 character id set by WebR3DAlign.  
% If Query.LoadFinal = 1, previous results are used; if 0, re-processed.
% Query.PostD is the discrepancy cutoff for post-processing
% Query.AnchorLength is the length of consecutive nt's required to be an
%    "anchor".  6 is the typical amount
% Query.SeedName is optional and is used to add an extension to the output
%    filenames
% Query.OutputFiles can be 0 or 1.  If 0, no output files are produced
% If Query.SaveEmptyAlignment = 1, alignments that do not find any
% correspondences are saved; if 0, not saved.
% seed1 and seed2 are optional seed alignments.

%Example 1: 
% Query.LoadFinal = 1;
% Query.Type='local';
% Query.PostD=3;
% Query.AnchorLength=6;
% Query.OutputFiles=0;
% Query.SaveEmptyAlignment = 1;
% R3DAlign('1j5e',{'A'},{'all'},'2avy',{'A'},{'all'},.5,8,10,'greedy',Query,Al1,Al2)

%Example 2: 
% Query.LoadFinal = 1;
% Query.Type='local';
% Query.PostD=3;
% Query.AnchorLength=6;
% Query.OutputFiles=1;
% Query.SaveEmptyAlignment = 1;
% [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign('4W21',{'8','5'},{'all','all'},'2QBG',{'B'},{'all'},{.4,.5,.5},{1,3,9},{200,70,20},{'greedy','greedy','greedy'},Query); 


%Differences between v2 and v1
% 1) Anchor code
% 2) Faster conversion of seed from cell
% 3) uses new rAugmentAlignment.m
% 4) uses rMakeEdgeMatrix.m

function [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,seed1,seed2)
try
if ~isfield(Query,'currIter')
   if nargin < 11
      Query.Type = 'local';
   end
   if ~isfield(Query,'Type')
      Query.Type = 'local';
   end
   if ~isfield(Query,'LoadFinal')
      Query.LoadFinal = 1;
   end
   if ~isfield(Query,'SeedName')
      Query.SeedName = '';
   end
   if ~isfield(Query,'AnchorLength')
      Query.AnchorLength = 15;
   end
   if ~isfield(Query,'SeqAlOutputFiles')
      Query.SeqAlOutputFiles = 0;
   end
   if ~isfield(Query,'OutputFiles')
      Query.OutputFiles = 1;
   end
   if ~isfield(Query,'ErrorMsg')
      Query.ErrorMsg = '';
   end
   if ~isfield(Query,'SaveEmptyAlignment')
      Query.SaveEmptyAlignment = 0;
   end
   if ~isfield(Query,'PostD')
      Query.PostD = 100;
   end   
   
%    addpath(genpath([pwd filesep 'FR3D']));
%Following was commented to be able to direct to specific R3DAlign folder
%manually
% addpath([pwd filesep 'R3DAlign']);  
   if ~(exist([pwd filesep 'PDBFiles']) == 7),        %#ok<*EXIST> % if directory doesn't yet exist
      mkdir([pwd filesep 'PDBFiles']);
   end
%    addpath([pwd filesep 'PDBFiles']);
   if ~(exist([pwd filesep 'PrecomputedData']) == 7),        % if directory doesn't yet exist
      mkdir([pwd filesep 'PrecomputedData']);
   end
%    addpath([pwd filesep 'PrecomputedData']);
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
      numNeigh{Query.currIter} = input('Enter the number of neighborhoods for each nucleotide: ');
   end
   if nargin < 7
      bandwidth{Query.currIter} = input('Enter the bandwidth for the seed alignment: ');
   end
   if nargin < 8
      cliqueMethod{Query.currIter} = input('Enter final clique method (Full or Greedy): ','s');
   end

   ErrorMsg='';
   AlignedNTs1=cell(2);
   AlignedNTs2=cell(2);

   NTList{1}=NTList1;
   NTList{2}=NTList2;
      
   % Check early (before pdb's are loaded) to see if alignment file already exists
   if ischar(File1) && ischar(File2)
      if isequal(File1,'uploaded')
         F1 = Query.UploadName1;
      else
         F1 = upper(File1);
      end
      if isequal(File2,'uploaded')
         F2 = Query.UploadName2;
      else
         F2 = upper(File2);
      end
      
      OutFilename = [F1 '(' Chain1{1}];
      for i=2:length(Chain1)
         OutFilename = [OutFilename '-' Chain1{i}]; %#ok<AGROW>  
      end
      OutFilename = [OutFilename ')'];
      
      for i=1:length(NTList1)
         if ~isequal(NTList1{i},'all')
            OutFilename = [OutFilename '(' NTList1{i} ')']; %#ok<AGROW>
         end
      end
      OutFilename = [OutFilename '--' F2 '(' Chain2{1}];
      for i=2:length(Chain2)
         OutFilename = [OutFilename '-' Chain2{i}]; %#ok<AGROW>  
      end
      OutFilename = [OutFilename ')'];
      for i=1:length(NTList2)
         if ~isequal(NTList2{i},'all')
            OutFilename = [OutFilename '(' NTList2{i} ')']; %#ok<AGROW>
         end
      end
      OutFilename = strrep(OutFilename, ':', '-');
      for i=1:length(discCut)
         OutFilename = [OutFilename '_d' num2str(discCut{i}) '_p' num2str(numNeigh{i}) '_B' num2str(bandwidth{i})]; %#ok<AGROW>
      end
      if ~isempty(Query.SeedName)
         OutFilename = [OutFilename '_' Query.SeedName];
      end
      OutFilename = strrep(OutFilename, '.', '');
      if exist(fullfile(pwd, 'R3D Align Output','Final Mat Files',F1, [OutFilename '.mat']))==2 && Query.LoadFinal==1
         disp('loading Final Alignment')
         load(fullfile(pwd, 'R3D Align Output', 'Final Mat Files',F1, [OutFilename '.mat']));
         return;
      else % Try if switching order of structures finds a file
         OutFilename = [F2 '('];
         for i=2:length(Chain2)
            OutFilename = [OutFilename '-' Chain2{i}]; %#ok<AGROW>  
         end
            OutFilename = [OutFilename ')'];
         for i=1:length(NTList2)
            if ~isequal(NTList2{i},'all')
               OutFilename = [OutFilename '(' NTList2{i} ')']; %#ok<AGROW>
            end
         end
         OutFilename = [OutFilename '--' F1 '('];
         for i=2:length(Chain1)
            OutFilename = [OutFilename '-' Chain1{i}]; %#ok<AGROW>  
         end
         OutFilename = [OutFilename ')'];
         for i=1:length(NTList1)
            if ~isequal(NTList1{i},'all')
               OutFilename = [OutFilename '(' NTList1{i} ')']; %#ok<AGROW>
            end
         end
         OutFilename = strrep(OutFilename, ':', '-');
         for i=1:length(discCut)
            OutFilename = [OutFilename '_d' num2str(discCut{i}) '_p' num2str(numNeigh{i}) '_B' num2str(bandwidth{i})]; %#ok<AGROW>
         end
         if ~isempty(Query.SeedName)
            OutFilename = [OutFilename '_' Query.SeedName];
         end
         OutFilename = strrep(OutFilename, '.', '');
         if exist(fullfile(pwd, 'R3D Align Output','Final Mat Files',F2, [OutFilename '.mat']))==2 && Query.LoadFinal==1
            disp('loading Final Alignment after switching structures')
            load(fullfile(pwd, 'R3D Align Output', 'Final Mat Files',F2, [OutFilename '.mat']));
            return;
         end
      end
   end

   try
      if ischar(File1),
        fprintf('Loading PDB Info...\n');
        if isequal(File1,'uploaded')
           Filename1 = Query.UploadName1;
           File1 = zAddNTData(Filename1,0);
           File1.Filename=Query.UploadName1;
        else
           Filename1 = upper(File1);
           File1 = zAddNTData(Filename1,0);
        end
      else %previously processed file was provided
         Filename1 = upper(File1.Filename);
      end
      Query.Filename1=Filename1;
   catch %#ok<CTCH>
      Query.ErrorMsg='Unable to load PDB file for Molecule 1.';
      [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
      return;
   end

   if File1.NumNT < 1
      Query.ErrorMsg='Unable to read PDB file for Molecule 1.';
      [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
      return;
   end

   try
      if ischar(File2),
         if isequal(File2,'uploaded')
            Filename2 = Query.UploadName2;
            File2 = zAddNTData(Filename2,0);
            File2.Filename=Query.UploadName2;
         else
            Filename2 = upper(File2);
            File2 = zAddNTData(Filename2,0);
         end
      else %previously processed file was provided
         Filename2 = upper(File2.Filename);
      end
      Query.Filename2=Filename2;
   catch %#ok<CTCH>
      Query.ErrorMsg='Unable to load PDB file for Molecule 2.';
      [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
      return;
   end

   if File2.NumNT < 1
      Query.ErrorMsg='Unable to read PDB file for Molecule 2.';
      [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
      return;
   end

   Indices1=[];
   Indices2=[];
%This code is all here to get the indices and determine if structures
% should be reversed before calling Iterative Align. 
   for i=1:length(NTList1)
      if isequal(Chain1{i},'all')
         Indices1=1:File1.NumNT;
      elseif isequal(NTList1{i},'all')
         chain=Chain1{i};
         c=cat(2,File1.NT.Chain);
      %Only use Model 1
         Mod=cat(2,File1.NT.ModelNum);
         tmpIndices = find(lower(c)==lower(chain) & Mod ==1);
         Indices1 = [Indices1 tmpIndices]; %#ok<AGROW>
      else
         if ~isempty(NTList1{i})
            tmpIndices = zIndexLookup(File1,NTList1{i},Chain1{i});
         else
            Query.ErrorMsg='No nucleotide numbers specified for Structure 1.';
            [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
            return 
         end
         Indices1 = [Indices1 tmpIndices]; %#ok<AGROW>
         if isempty(Indices1)
            Query.ErrorMsg='Invalid nucleotide numbers specified for Structure 1.';
            [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
            return; 
         end
      end
   end

   for i=1:length(NTList2)
      if isequal(Chain2{i},'all')
         Indices2=1:File2.NumNT;
      elseif isequal(NTList2{i},'all')
         chain=Chain2{i};
         c=cat(2,File2.NT.Chain);
         Mod=cat(2,File2.NT.ModelNum);
         tmpIndices = find(lower(c)==lower(chain) & Mod ==1);
         Indices2 = [Indices2 tmpIndices]; %#ok<AGROW>
      else
         if ~isempty(NTList2{i})
            tmpIndices = zIndexLookup(File2,NTList2{i},Chain2{i});
         else
            Query.ErrorMsg='No nucleotide numbers specified for Structure 2.';
            [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
            return;
         end
         Indices2 = [Indices2 tmpIndices]; %#ok<AGROW>
         if isempty(Indices2)
            Query.ErrorMsg='Invalid nucleotide numbers specified for Structure 2.';
            [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
            return; 
         end
      end
   end
   
   if length(Indices1)<4     
      Query.ErrorMsg='Structure 1 must contain at least four nucleotides.';
      [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
      return;
   elseif length(Indices2)<4     
      Query.ErrorMsg='Structure 2 must contain at least four nucleotides.';
      [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
      return;
   end

%7/8/2014 Smaller structure always first - faster and better seed alignment
%augmentation.  Example of benefit: 1n36A vs. 1sloA
   if length(Indices1) > length(Indices2)
      temp = File2;
      File2 = File1;
      File1 = temp;
      temp=Query.Filename2;
      Query.Filename2=Query.Filename1;
      Query.Filename1=temp;
      temp = Chain2;
      Chain2 = Chain1;
      Chain1 = temp;
      temp = NTList2;
      NTList2 = NTList1;
      NTList1 = temp;
      temp = Indices2;
      Indices2 = Indices1;
      Indices1 = temp;
%       if exist('seed1') || ~isempty(seed1)
      if exist('seed1')
         temp = seed2;
         seed2 = seed1;
         seed1 = temp;
      end
   end

   if iscell(discCut)
      if ~exist('seed1') || isempty(seed1)
         [AlignedNTs1,AlignedNTs2,ErrorMsg] = IterativeAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query);
      else
         [AlignedNTs1,AlignedNTs2,ErrorMsg] = IterativeAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,seed1,seed2);
      end
      return;
   else  %there will be only one iteration and convert variables into cell format
      Query.currIter = 1; %Iteration Number
      discCut={discCut};
      numNeigh={numNeigh};
      bandwidth={bandwidth};
      cliqueMethod={cliqueMethod};
   end
end

Indices1=[];
Indices2=[];
ErrorMsg='';
AlignedNTs1=cell(2);
AlignedNTs2=cell(2);
NTList{1}=NTList1;
NTList{2}=NTList2;
OutFilename = [Query.Filename1 '(' Chain1{1}];
NeighAFilename = ['Neighborhoods' filesep Query.Filename1 '(' Chain1{1}];

for i=2:length(Chain1)
   OutFilename = [OutFilename '-' Chain1{i}]; %#ok<AGROW>  
   NeighAFilename = [NeighAFilename '-' Chain1{i}]; %#ok<AGROW> 
end
OutFilename = [OutFilename ')'];  
NeighAFilename = [NeighAFilename ')']; 

for i=1:length(NTList1)
   if ~isequal(NTList1{i},'all')
      OutFilename = [OutFilename '(' NTList1{i} ')']; %#ok<AGROW>
      NeighAFilename = [NeighAFilename '(' NTList1{i} ')']; %#ok<AGROW>
   end
   if isequal(Chain1{i},'all')
      Indices1=1:File1.NumNT;
   elseif isequal(NTList1{i},'all')
      chain=Chain1{i};
      c=cat(2,File1.NT.Chain);
   %Only use Model 1
      Mod=cat(2,File1.NT.ModelNum);
      tmpIndices = find(lower(c)==lower(chain) & Mod ==1);
      Indices1 = [Indices1 tmpIndices]; %#ok<AGROW>
   else
      if ~isempty(NTList1{i})
         tmpIndices = zIndexLookup(File1,NTList1{i},Chain1{i});
      else
         Query.ErrorMsg='No nucleotide numbers specified for Structure 1.';
         [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
         return 
      end
      Indices1 = [Indices1 tmpIndices]; %#ok<AGROW>
      if isempty(Indices1)
        Query.ErrorMsg='Invalid nucleotide numbers specified for Structure 1.';
        [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
        return; 
      end
   end
end
NeighAFilename = [NeighAFilename '_' num2str(numNeigh{Query.currIter}) '.mat'];
NeighAFilename = strrep(NeighAFilename, ':', '-');
NeighAFilename = [pwd filesep NeighAFilename];

OutFilename = [OutFilename '--' Query.Filename2 '(' Chain2{1}];
NeighBFilename = ['Neighborhoods' filesep Query.Filename2 '(' Chain2{1}];


for i=2:length(Chain2)
   OutFilename = [OutFilename '-' Chain2{i}]; %#ok<AGROW>  
   NeighBFilename = [NeighBFilename '-' Chain2{i}]; %#ok<AGROW> 
end
OutFilename = [OutFilename ')'];  
NeighBFilename = [NeighBFilename ')'];

for i=1:length(NTList2)
   if ~isequal(NTList2{i},'all')
      OutFilename = [OutFilename NTList2{i}]; %#ok<AGROW>
      NeighBFilename = [NeighBFilename NTList2{i}]; %#ok<AGROW>
   end
   if isequal(Chain2{i},'all')
      Indices2=1:File2.NumNT;
   elseif isequal(NTList2{i},'all')
      chain=Chain2{i};
      c=cat(2,File2.NT.Chain);
      Mod=cat(2,File2.NT.ModelNum);
      tmpIndices = find(lower(c)==lower(chain) & Mod ==1);
      Indices2 = [Indices2 tmpIndices]; %#ok<AGROW>
   else
      if ~isempty(NTList2{i})
         tmpIndices = zIndexLookup(File2,NTList2{i},Chain2{i});
      else
         Query.ErrorMsg='No nucleotide numbers specified for Structure 2.';
         [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
         return 
      end
      Indices2 = [Indices2 tmpIndices]; %#ok<AGROW>
      if isempty(Indices2)
        Query.ErrorMsg='Invalid nucleotide numbers specified for Structure 2.';
        [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
        return; 
      end
   end
end
if length(Indices1)<4     
   Query.ErrorMsg='Structure 1 must contain at least four nucleotides.';
   [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
   return;
elseif length(Indices2)<4     
   Query.ErrorMsg='Structure 2 must contain at least four nucleotides.';
   [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
   return;
end
OutFilename = strrep(OutFilename, ':', '-');
ShortOutFilename=OutFilename;
for i=1:Query.currIter
   OutFilename = [OutFilename '_d' num2str(discCut{i}) '_p' num2str(numNeigh{i}) '_B' num2str(bandwidth{i})]; %#ok<AGROW>
end
if ~isempty(Query.SeedName)
  OutFilename = [OutFilename '_' Query.SeedName];
end
OutFilename = strrep(OutFilename, '.', '');

NeighBFilename = [NeighBFilename '_' num2str(numNeigh{Query.currIter}) '.mat'];
NeighBFilename = strrep(NeighBFilename, ':', '-');
NeighBFilename = [pwd filesep NeighBFilename];

%7/8/2014 Smaller structure always first - faster and better seed alignment
%augmentation.  Example of benefit: 1n36A vs. 1sloA

% if length(Indices1) > length(Indices2)
%    temp = NeighBFilename;
%    NeighBFilename = NeighAFilename;
%    NeighAFilename = temp;
%    temp = File2;
%    File2 = File1;
%    File1 = temp;
%    temp = Chain2;
%    Chain2 = Chain1;
%    Chain1 = temp;
%    temp = NTList2;
%    NTList2 = NTList1;
%    NTList1 = temp;
%    temp = Indices2;
%    Indices2 = Indices1;
%    Indices1 = temp;
%    if exist('seed1') || ~isempty(seed1)
%       temp = seed2;
%       seed2 = seed1;
%       seed1 = temp;
%    end
% end
   
if exist(fullfile(pwd, 'R3D Align Output','Final Mat Files', Query.Filename1, [OutFilename '.mat']))==2 && Query.LoadFinal==1 
   disp('loading Final Alignment')
   load(fullfile(pwd, 'R3D Align Output', 'Final Mat Files', Query.Filename1, [OutFilename '.mat']));
else
   disp('not loading Final Alignment')
   if numNeigh{Query.currIter} > nchoosek(length(Indices1)-1,3)     
      Query.ErrorMsg='Neighborhood parameter too large based on size of structure 1.';
      [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
      return;
   elseif numNeigh{Query.currIter} > nchoosek(length(Indices2)-1,3)
      Query.ErrorMsg='Neighborhood parameter too large based on size of structure 2.';
      [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
      return;
   end

   disp('get neighborhoods')

   if exist(NeighAFilename)==2
      disp(['loading ' NeighAFilename])
      load(NeighAFilename);
      AQuads=Quads; %#ok<NODEF>
      clear Quads;
      c = cat(1,File1.NT(Indices1).Center);           
      try
         File1.Distance = full(zMutualDistance(c,Inf));
      catch %#ok<CTCH>
         Query.ErrorMsg='Error in zMutualDistance.';
         [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
         return;
      end
      A = triu(File1.Distance); % will need matrix A later
   else
      disp(['not loading ' NeighAFilename])
      maxdist=15;
      getMoreNeigh = true;
      while getMoreNeigh == true;
         c = cat(1,File1.NT(Indices1).Center);           % nucleotide centers
         File1.Distance = full(zMutualDistance(c,Inf));  % Distance matrix for A
         A = triu(File1.Distance);                       % Distance matrix is symmetrical;
         [rA, cA] = find(A<maxdist & A>0);
         [S1, IX] = sort(rA);
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
      end
      AQuads = rGetQuads(AQuads,A,numNeigh{Query.currIter});
      Quads=AQuads; 
      save(NeighAFilename, 'Quads')
      clear Quads;
   end

   if exist(NeighBFilename)==2
      disp(['loading ' NeighBFilename])
      load(NeighBFilename);
      BQuads=Quads;
      clear Quads;
      d = cat(1,File2.NT(Indices2).Center);           % will need matrix B later
      File2.Distance = full(zMutualDistance(d,Inf));
      B = triu(File2.Distance);
   else
      disp(['not loading ' NeighBFilename])
      maxdist=15;
      getMoreNeigh = true;
      while getMoreNeigh == true;   
         d = cat(1,File2.NT(Indices2).Center);           % nucleotide centers
         File2.Distance = full(zMutualDistance(d,Inf));  % Distance matrix for B
         B = triu(File2.Distance);                       % don't want duplicate entries returned in next lines
         [rB, cB] = find(B<maxdist & B>0);
         [S2, IX] = sort(rB);
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

   %Web server does not allow for seed uploading at this time so following
   %if block should not be utilized ever.
   if isequal(Query.Type,'web') && ~isequal(Query.SeedName,'')
      FASTA = zReadFASTA(['seed_upload_file' '.fasta']);
      delete(['./SeedAlignments/' Query.Name '.fasta']);
      File1Sequence=cat(2,File1.NT(Indices1).Base);
      File2Sequence=cat(2,File2.NT(Indices2).Base);

      if abs(length(FASTA(1).Sequence)-length(File1Sequence))>0 || abs(length(FASTA(2).Sequence)-length(File2Sequence))>0
         Query.ErrorMsg='Seed alignment does not correspond with pdb files.';
         [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
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
   elseif ~exist('seed1') || isempty(seed1)
      if exist(fullfile(pwd, 'Sequence Alignments', [ShortOutFilename '.mat']))==2
         disp('loading Sequence Alignment')
         load(fullfile(pwd, 'Sequence Alignments', [ShortOutFilename '.mat']));
      else
         disp('not loading Sequence Alignment')
         [align1, align2, charAlign1, charAlign2] = rGapNW(File1,Indices1,File2,Indices2,.999,7,0.25); %#ok<NASGU,ASGLU>
         save(fullfile(pwd, 'Sequence Alignments', [ShortOutFilename '.mat']), 'align1', 'align2', 'charAlign1', 'charAlign2')
      end
   else
%       for i=1:length(seed1)
%          if iscell(seed1)
%             tmp=zIndexLookup(File1,seed1{i,1},seed1{i,2},0);
%             align1(i)=find(Indices1==tmp); %#ok<AGROW>
%             tmp=zIndexLookup(File2,seed2{i,1},seed2{i,2},0);
%             align2(i)=find(Indices2==tmp); %#ok<AGROW>
%          elseif isvector(seed1)
%             align1=seed1;
%             align2=seed2;
%          end
%       end

%This code should be faster:
      if iscell(seed1) 
         tmp1=zIndexLookup(File1,seed1(:,1),seed1(:,2),0);
         [tf, align1] = ismember(tmp1,Indices1);
         tmp2=zIndexLookup(File2,seed2(:,1),seed2(:,2),0);
         [tf, align2] = ismember(tmp2,Indices2);
      elseif isvector(seed1)
         align1=seed1;
         align2=seed2;
      end  
         
   end

%Files for the sequence alignment are created and moved to the
%'Sequence Alignments' folder.
   if ~isequal(Query.Type,'web') && exist(fullfile(pwd, 'Sequence Alignments',ShortOutFilename)) ~= 7
      I1=Indices1(align1);
      I2=Indices2(align2);
      
      clf
      if Query.SeqAlOutputFiles == 1
         rAlignmentSpreadsheet(File1,Indices1,File2,Indices2,I1,I2,ShortOutFilename,ErrorMsg);
         try
            [AAA,BBB] = rBarDiagram(File1,Indices1,File2,Indices2,I1,I2,ShortOutFilename,'R3D Align');
         catch %#ok<CTCH>
             fprintf('Bar diagram could not be generated for sequence alignment\n');
         end
      
         View = [1 1 1 1 0 0 0];
         try
            m1 = zBarDiagramInteractions(File1,Indices1,AAA,View,'above');
            m2 = zBarDiagramInteractions(File2,Indices2,BBB,View,'below');
            if m1(4) ~= 0 && m2(3) ~= 0
               axis([0 20 m2(3) m1(4)])
            end
            saveas(gcf,[ShortOutFilename '_int'],'pdf')
         catch %#ok<CTCH>
            fprintf('Interaction bar diagram could not be generated\n');
         end
         rWriteAlignmentFasta(File1,Indices1,File2,Indices2,I1,I2,NTList,ShortOutFilename);
         rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,I1,I2,NTList,ShortOutFilename);
         movefile([pwd filesep ShortOutFilename '*'], fullfile(pwd, 'Sequence Alignments', ShortOutFilename));
      end
      clear I1 I2;
   end

%*************************************
   if length(align1)>4
      D=rFindAlignmentDiscrepancies(File1,Indices1(align1),File2,Indices2(align2),'nearest4');
      D=D';
   else
      D=10*ones(1,length(align1));
   end
%*************************************
 
 % Augment alignment with near matches
% For each nt in A aligned with a gap, the nearest aligned nt in B is found
% to use as the center of the band for that nt
%*************************************
   if length(Indices1) ~= length(align1)  % if not all nts had a correspondence in seq alignment
      [align1,align2,D] = rAugmentAlignment(Indices1,Indices2,align1,align2,D); 
   end
  
   D=D<.5;
   R=D*0;
   Anchor=D*0+1;
   R(1)=D(1);
   for i=2:length(D)
      R(i)=(R(i-1)+1)*D(i);  %R(i) is the number of consecutive anchor points leading to i
   end
   currMax = max(R);
   while currMax > Query.AnchorLength
       loc=find(R==currMax);
       Anchor(loc-(currMax-2):loc-1) = 0;  %Assign anchor designations, but not to end boundary nt's
       R(loc-(currMax-1):loc) = 0;
       currMax = max(R);
   end
   
%  Anchor = Anchor == 0;    %Use 0 for anchored and 1 for non-anchored
   LB=Anchor*0;
   RB=Anchor*0;
   LBtmp=0;

   for i=1:length(Anchor)
      if Anchor(i)== 1
         Next=find(Anchor(i+1:end)==0,1,'first');
         LB(i)=LBtmp+1;
         if ~isempty(Next)
            RB(i)=align2(i + Next)-1;
         else
            RB(i)=length(Indices2);
         end
      else
         LBtmp=align2(i);
         LB(i)=align2(i);
         RB(i)=align2(i);
      end
   end
   disp(['Number of Anchor Points:' num2str(sum(Anchor==0)) '/' num2str(length(Anchor))]);
  
   SA=round(bandwidth{Query.currIter}/2);
   LBtmp=align2-SA;
   RBtmp=align2+SA;
   LB=max(LB,LBtmp);
   RB=min(RB,RBtmp);
%*************************************
%    SB=[];
   VMI=[];         %This will be a matrix containing the vertices that are created
                  %eg. if VMI(1) =[2 3], there was a vertex created for the
                  %2nd A neighborhood and 3rd B neighborhood
   nA=length(AQuads(:,1));
%    nB=length(BQuads(:,1));

   List=cell(length(A),length(B));     %List{i,j} is a vector indicating the vertices that align
                                    %nucleotide i of A with nucleotide j of B
   Len=zeros(length(A),length(B),'uint16');

% Cutoff values for screening criterion
   cut2=4*4*2*discCut{Query.currIter}*discCut{Query.currIter};
   cut3=4*4*3*discCut{Query.currIter}*discCut{Query.currIter};
   cut4=4*4*4*discCut{Query.currIter}*discCut{Query.currIter};
%    SA=round(bandwidth{Query.currIter}/2);
%    fprintf('Making Local Alignment Graph...\n');
% Determine locations of vertices

   for i = 1:nA
     if Anchor(AQuads(i,1))==1 || Anchor(AQuads(i,2))==1 || Anchor(AQuads(i,3))==1 || Anchor(AQuads(i,4))==1
%6/9/2014 Test this sometime for comparison: if sum(Anchor(AQuads(i,:))>0

%       SB(1)=align2(AQuads(i,1));
%       SB(2)=align2(AQuads(i,2));
%       SB(3)=align2(AQuads(i,3));
%       SB(4)=align2(AQuads(i,4));
 
T1=find(ismember(BQuads(:,1),LB(AQuads(i,1)):RB(AQuads(i,1))));
T2=find(ismember(BQuads(:,2),LB(AQuads(i,2)):RB(AQuads(i,2))));
T3=find(ismember(BQuads(:,3),LB(AQuads(i,3)):RB(AQuads(i,3))));
T4=find(ismember(BQuads(:,4),LB(AQuads(i,4)):RB(AQuads(i,4))));
T=intersect(intersect(T1,T2),intersect(T3,T4));
%       for j = 1:nB
     for j=T'
         if isempty(j)
             break;
         end
%          if BQuads(j,1) >= LB(AQuads(i,1)) && BQuads(j,1) <= RB(AQuads(i,1)) && ...
%             BQuads(j,2) >= LB(AQuads(i,2)) && BQuads(j,2) <= RB(AQuads(i,2)) && ...     
%             BQuads(j,3) >= LB(AQuads(i,3)) && BQuads(j,3) <= RB(AQuads(i,3)) && ...     
%             BQuads(j,4) >= LB(AQuads(i,4)) && BQuads(j,4) <= RB(AQuads(i,4))     
%          if abs(SB(1)-BQuads(j,1))<SA && abs(SB(2)-BQuads(j,2))<SA  && abs(SB(3)-BQuads(j,3))<SA && abs(SB(4)-BQuads(j,4))<SA
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
                  Len(AQuads(i,1),BQuads(j,1))=length(List{AQuads(i,1),BQuads(j,1)});
                  Len(AQuads(i,2),BQuads(j,2))=length(List{AQuads(i,2),BQuads(j,2)});
                  Len(AQuads(i,3),BQuads(j,3))=length(List{AQuads(i,3),BQuads(j,3)});
                  Len(AQuads(i,4),BQuads(j,4))=length(List{AQuads(i,4),BQuads(j,4)});
               end
            end
%          end
      end
     end
   end

   if isempty(VMI)
      AlignedIndices1 = [];
      AlignedIndices2 = [];
      AlignedNTs1=[];
      AlignedNTs2=[];
      if Query.OutputFiles == 1
         WriteOutput(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename,ShortOutFilename,ErrorMsg,Query)
      end
      if Query.SaveEmptyAlignment == 1
         save(fullfile(pwd, 'R3D Align Output', 'Final Mat Files', Query.Filename1, [OutFilename '.mat']),'Indices1','Indices2','AlignedIndices1','AlignedIndices2','AlignedNTs1','AlignedNTs2');
      end
      return;
   end
   
try
   [EdgeMat]=rMakeEdgeMatrix(VMI,List,Len);
% % % catch %#ok<CTCH>
% % %   Query.ErrorMsg = ['Edge Matrix Memory Error in Iteration ' num2str(Query.currIter) '.'];
% % %   [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
% % %   return;
% % % end  
catch ME
       for k=1:length(ME.stack)
          ME.stack(k)
       end
       Query.ErrorMsg = ['Edge Matrix Memory Error in Iteration ' num2str(Query.currIter) '.'];
       [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
       return;
end
   clear List;
   clear A;
   clear B;
%*************************************
Diff=[Anchor(1) diff(Anchor)];
Starts=find(Diff==1);
Ends=find(Diff==-1)-1;
if Anchor(end) == 1
   Ends(end+1)=length(Anchor);
end
Segments = [Starts' Ends'];
Indnums1 = align1(Anchor==0)';
Indnums2 = align2(Anchor==0)';

  AlignedIndices1 = Indices1(Indnums1);
  AlignedIndices2 = Indices2(Indnums2);
%     [AAA,BBB] = rBarDiagram(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,[OutFilename 'Anchor' num2str(Query.currIter)],'R3D Align');
%*******************************
% fprintf('Here we go...\n');
for i = 1:length(Segments(:,1))
    [t1, t2]= find(ismember(AQuads(VMI(:,1),:), Starts(i):Ends(i))); %#ok<NASGU>
    t1=sort(unique(t1));
    EM=EdgeMat(t1,t1);
    V=VMI(t1,:);

   maximalClique = rGetMaximalClique(EM);
   CLB = length(maximalClique);  %Current lower bound

   refineMore=true;
   numLoop=0;
%    fprintf('Preprocessing...\n');
   while refineMore==true
      numLoop=numLoop+1;
      D=[];
      C=colorVertices(EM);       %vector of length=length(EM(:,1)) with color class for each vertex
      if ~isempty(EM)
         for p = 1:length(EM(:,1))
            if ~any(maximalClique==p)
               CX = C(EM(p,:)==0);               %Colors of neighbors
               if length(unique(CX))<CLB-1
                  D=[D p]; %#ok<AGROW>
               end
            end
         end
      
         if length(D)/length(EM(:,1)) < .2 || numLoop == 2
            refineMore=false;
         end
      else
         refineMore=false;
      end
      V(D,:)=[]; 
      EM(D,:)=[];
      EM(:,D)=[];
   end
%    fprintf('Finding Clique...\n');

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


% Convert clique into corresponding alignment
   Z=V(MC,:);
   tmpAlmnt = [AQuads(Z(:,1),1) BQuads(Z(:,2),1);
               AQuads(Z(:,1),2) BQuads(Z(:,2),2);
               AQuads(Z(:,1),3) BQuads(Z(:,2),3);
               AQuads(Z(:,1),4) BQuads(Z(:,2),4)];  %unsorted, with duplicates
   [S1, IX] = sort(tmpAlmnt(:,1));
   tmpAlmnt(:,1)=S1;
   tmpAlmnt(:,2) = tmpAlmnt(IX,2);    %sorted
   [b, m] = unique(tmpAlmnt(:,1));
   Indnums1tmp=b;
   Indnums2tmp=tmpAlmnt(m,2);            %no duplicates
   T=~ismember(Indnums1tmp,Starts(i):Ends(i));
   Indnums1tmp(T)=[];
   Indnums2tmp(T)=[];
   Indnums1=[Indnums1; Indnums1tmp];
   Indnums2=[Indnums2; Indnums2tmp];
end



%The following are the indices of aligned nucleotides

   AlignedIndices1 = Indices1(Indnums1);
   AlignedIndices2 = Indices2(Indnums2);
   OrigIndices1 = AlignedIndices1;
   c=cat(2,File1.NT(AlignedIndices1).Chain);
   for i=1:length(Chain1)
      I = find(lower(c)==lower(Chain1{i}));
      AlignedIndices1(I) = AlignedIndices1(I)+100000*(i-1);
   end
   [~, ind] = sort(AlignedIndices1);
   AlignedIndices1 = OrigIndices1(ind);
   AlignedIndices2 = AlignedIndices2(ind);
  
   
if length(Indnums1)>4
%       sum(Anchor)
      D=rFindAlignmentDiscrepancies(File1,AlignedIndices1,File2,AlignedIndices2,'nearest4');
      disp(['Number removed via post-processing: ' num2str(sum(D>Query.PostD)) '/' num2str(length(D))]);
      AlignedIndices1(D>Query.PostD)=[];
      AlignedIndices2(D>Query.PostD)=[];
end
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
   if ~(exist(fullfile(pwd, 'R3D Align Output', 'Final Mat Files', Query.Filename1)) == 7),        % if directory doesn't yet exist
      mkdir(fullfile(pwd, 'R3D Align Output', 'Final Mat Files', Query.Filename1));
   end
   save(fullfile(pwd, 'R3D Align Output', 'Final Mat Files', Query.Filename1, [OutFilename '.mat']),'Indices1','Indices2','AlignedIndices1','AlignedIndices2','AlignedNTs1','AlignedNTs2');
end

if Query.OutputFiles == 1
   fprintf('Producing Alignment Output...\n');
   WriteOutput(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename,ShortOutFilename,ErrorMsg,Query)
end
% fprintf('Done\n');
catch ME
  ME.identifier
  ME.message
  if (strcmp(ME.identifier,'MATLAB:nomem'))
     Query.ErrorMsg = ['Memory Error in Iteration ' num2str(Query.currIter) '.'];
     [AlignedNTs1, AlignedNTs2, ErrorMsg] = HandleError(Query);
     return;
  else
     for k=1:length(ME.stack)
       ME.stack(k)
     end
  end
  ErrorMsg = ME.message;
  return;
end
end

function [AlignedNTs1,AlignedNTs2,ErrorMsg] = HandleError(Query)
   AlignedNTs1{1,1}=[];
   AlignedNTs2{1,1}=[];
   ErrorMsg=Query.ErrorMsg;
   if strcmpi(Query.Type,'local')
      fprintf('Error: %s \n',ErrorMsg);
   elseif isequal(Query.Type,'web')
      save([pwd filesep Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
   end
end

function WriteOutput(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename,ShortOutFilename,ErrorMsg,Query)
if isequal(Query.Type,'web')
   if strcmp(ErrorMsg,'')
      clf
      [AAA,BBB] = rBarDiagram(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,Query.Name,'R3D Align');
      View = [1 1 1 1 0 0 0];
      m1 = zBarDiagramInteractions(File1,Indices1,AAA,View,'above');
      m2 = zBarDiagramInteractions(File2,Indices2,BBB,View,'below');
      if m1(4) ~= 0 && m2(3) ~= 0
         axis([0 20 m2(3) m1(4)])
      end
      saveas(gcf,[Query.Name '_int'],'pdf')
      saveas(gcf,[Query.Name '_int'],'png')      
      
      rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
      rWriteAlignmentFasta(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
      VP.Write=1;
      rSuperimposeNucleotides(File1,[AlignedIndices1 setdiff(Indices1,AlignedIndices1)],File2,[AlignedIndices2 setdiff(Indices2,AlignedIndices2)],VP,length(AlignedIndices1),Query.Name);
      FL = rAlignmentSpreadsheet(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,Query.Name,ErrorMsg);
      rWriteSummaryStatistics(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,Query.Name,FL);
      save([pwd filesep Query.Name '.mat'], 'AlignedIndices1', 'AlignedIndices2', 'ErrorMsg', 'Query');
   end
elseif exist(fullfile(pwd, 'R3D Align Output', OutFilename)) ~= 7 %if folder does not exist
   clf
   try
       [AAA,BBB] = rBarDiagram(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,OutFilename,'R3D Align');
   catch ME
       for k=1:length(ME.stack)
          ME.stack(k)
       end
       fprintf('Bar diagram could not be generated\n');
   end
   View = [1 1 1 1 0 0 0];
   try
       m1 = zBarDiagramInteractions(File1,Indices1,AAA,View,'above');
       m2 = zBarDiagramInteractions(File2,Indices2,BBB,View,'below');
       if m1(4) ~= 0 && m2(3) ~= 0
          axis([0 20 m2(3) m1(4)])
       end
       saveas(gcf,[OutFilename '_int'],'pdf')
   catch ME
       for k=1:length(ME.stack)
          ME.stack(k)
       end
       fprintf('Interaction bar diagram could not be generated\n');
   end
   rWriteAlignmentFasta(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename);
   FL=rAlignmentSpreadsheet(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,OutFilename,ErrorMsg);
   rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,OutFilename);
   movefile([pwd filesep OutFilename '*'], fullfile(pwd, 'R3D Align Output', OutFilename));
   try
       rAnalyzeAlignmentNew(File1,File2,Indices1,Indices2,AlignedIndices1,AlignedIndices2,OutFilename,ShortOutFilename,ErrorMsg,Query); 
   catch ME
       for k=1:length(ME.stack)
          ME.stack(k)
       end
       fprintf('Alignment could not be analyzed\n');
   end
     
   if exist('Query', 'var') && isfield(Query, 'Name')
       VP.Write=1;
       rSuperimposeNucleotides(File1,[AlignedIndices1 setdiff(Indices1,AlignedIndices1)],File2,[AlignedIndices2 setdiff(Indices2,AlignedIndices2)],VP,length(AlignedIndices1),Query.Name);
   end
end                  
end

function [AlignedNTs1,AlignedNTs2,ErrorMsg] =IterativeAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,seed1,seed2)
   
   L=length(discCut); %number of iterations for the alignment
   Query.currIter = 1; %current Iteration
   if ~exist('seed1') || isempty(seed1)
      [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query);
   else
      [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,seed1,seed2);
   end
   ct=2;  %Just completed first iteration, now starting second
   while ct <= L
      if ~strcmpi(ErrorMsg,'') %there is an Error, so stop execution
         break;
      else
         Query.currIter = ct;
         Query.ErrorMsg=ErrorMsg;
         if ~isempty(AlignedNTs1) %must check to see if previous alignment aligned any nt's to prevent error of AlignedNTS1 being empty
            [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,AlignedNTs1,AlignedNTs2);
         else
            [AlignedNTs1,AlignedNTs2,ErrorMsg] = R3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query);
         end
         ct = ct + 1;
      end
   end
end
