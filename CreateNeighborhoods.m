% FileList is  a cell array of pdb ids
% Filelist = {'1s7s','2qbg','2j01'}

function [] = CreateNeighborhoods(FileList)

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
if ~(exist([pwd filesep 'Neighborhoods']) == 7),        % if directory doesn't yet exist
   mkdir([pwd filesep 'Neighborhoods']);
end

for i=1:length(FileList)
   Filename = upper(FileList{i});
   File = zAddNTData(Filename,0); 
   NeighFilename1 = fullfile(pwd, 'Neighborhoods', Filename);
   Chains = unique(cat(2,File.NT.Chain));
   for j=1:length(Chains)
      NeighFilename2 = [NeighFilename1 '(' Chains(j) ')' 'all']; %#ok<AGROW>
      c=cat(2,File.NT.Chain);
      Indices = find(lower(c)==lower(Chains(j)));
      for k=1:11
         NeighFilename3 = [NeighFilename2 '_' num2str(k) '.mat']; %#ok<AGROW>
         if k > nchoosek(length(Indices)-1,3)
            disp(['Too many neighborhoods for a small structure...' NeighFilename3]) 
         elseif exist(NeighFilename3)~=2
            maxdist=15;
            getMoreNeigh = true;
            while getMoreNeigh == true;
               disp(['Creating ' NeighFilename3])
               c = cat(1,File.NT(Indices).Center);           % nucleotide centers
               File.Distance = full(zMutualDistance(c,Inf));  % Distance matrix for A
               A = triu(File.Distance);                       % Distance matrix is symmetrical;
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
               AQuads = rGetQuads(AQuads,A,k);
               Quads=AQuads;  %#ok<NASGU>
               save(NeighFilename3, 'Quads')
               clear Quads;
            end
         else
            disp([NeighFilename3 ' already exists']) 
         end
      end
   end
end