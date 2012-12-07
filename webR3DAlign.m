function [AlignedNTs1,AlignedNTs2,ErrorMsg] = webR3DAlign(File1,Chain1,NTList1,File2,Chain2,NTList2,discCut,numNeigh,bandwidth,cliqueMethod,Query,seed1,seed2)


%     webroot = '/Servers/rna.bgsu.edu/WebR3DAlign';
%     fr3dpath='/Servers/rna.bgsu.edu/WebR3DAlign/FR3D';

    webroot = [filesep 'Servers' filesep 'rna.bgsu.edu' filesep 'r3dalign_dev'];
%     fr3dpath= [pwd filesep 'FR3D'];
%     addpath(genpath(fr3dpath));
    
% fr3dpath=pwd;
%     webroot=pwd;

path(path,[pwd filesep 'FR3D' filesep 'FR3DSource']);

if ~(exist('PDBFiles') == 7),        % if directory doesn't yet exist
  mkdir('PDBFiles');
end
path(path,[pwd filesep 'PDBFiles']);

if ~(exist('PrecomputedData') == 7),        % if directory doesn't yet exist
  mkdir('PrecomputedData');
end
path(path,[pwd filesep 'PrecomputedData']);

path(path,[pwd filesep 'R3DAlign']);

ResultFolder = fullfile(webroot, 'data', 'results', Query.Name);
if~(exist(ResultFolder,'dir')),
   mkdir(ResultFolder);
end
addpath(ResultFolder);




% try
ErrorMsg='';
AlignedNTs1=cell(2);
AlignedNTs2=cell(2);

NTList{1}=NTList1;
NTList{2}=NTList2;

if ischar(File1),
  Filename = File1;
end

try
if isequal(File1,'uploaded')
   Filename = [Query.Name '_1'];
   File1 = zGetNTData(Filename,0);
   File1.Filename=Query.UploadName1;
else
  Filename = File1;
  File1 = zGetNTData(Filename,0);
end
catch
  ErrorMsg='ERROR: Unable to load PDB file for Molecule 1.';
  save(fullfile(ResultFolder, [Query.Name '.mat']), 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
  return;
end

if isempty(File1.Backbone)
   ErrorMsg='ERROR: Unable to load PDB file for Molecule 1.';
   save(fullfile(ResultFolder, [Query.Name '.mat']), 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
   return;
end

if ischar(File2),
  Filename = File2;
end

try
if isequal(File2,'uploaded')
   Filename = [Query.Name '_2'];
   File2 = zGetNTData(Filename,0);
   File2.Filename=Query.UploadName2;
else
  Filename = File2;
  File2 = zGetNTData(Filename,0);
end
catch
  ErrorMsg='ERROR: Unable to load PDB file for Molecule 2.';
  save(fullfile(ResultFolder, [Query.Name '.mat']), 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
  return;
end

if isempty(File2.Backbone)
   ErrorMsg='ERROR: Unable to load PDB file for Molecule 2.';
   save(fullfile(ResultFolder, [Query.Name '.mat']), 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
   return;
end

Indices1=[];
for i=1:length(NTList1)
   if isequal(Chain1{i},'all')
      Indices1=1:File1.NumNT;
   elseif isequal(NTList1{i},'all')
      chain=Chain1{i};
      c=cat(2,File1.NT.Chain);
      tmpIndices = find(lower(c)==lower(chain));
      Indices1 = [Indices1 tmpIndices];
   else
      tmpIndices = zIndexLookup(File1,NTList1{i},Chain1{i});
      Indices1 = [Indices1 tmpIndices];
   end
end
Indices2=[];
for i=1:length(NTList2)
   if isequal(Chain2{i},'all')
      Indices2=1:File2.NumNT;
   elseif isequal(NTList2{i},'all')
      chain=Chain2{i};
      c=cat(2,File2.NT.Chain);
      tmpIndices = find(lower(c)==lower(chain));
      Indices2 = [Indices2 tmpIndices];
   else
      tmpIndices = zIndexLookup(File2,NTList2{i},Chain2{i});
      Indices2 = [Indices2 tmpIndices];
   end
end  
if numNeigh > nchoosek(length(Indices1)-1,3)
    ErrorMsg='ERROR: Neighborhood parameter too high based on size of structure.';
    save(fullfile(ResultFolder, [Query.Name '.mat']), 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
    return;
elseif numNeigh > nchoosek(length(Indices2)-1,3)
    ErrorMsg='ERROR: Neighborhood parameter too high based on size of structure.';
    save(fullfile(ResultFolder, [Query.Name '.mat']), 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
    return;
end
if ~isequal(Query.SeedName,'')
   FASTA = zReadFASTA([Query.Name '.fasta']);
   delete(['./SeedAlignments/' Query.Name '.fasta']);
   File1Sequence=cat(2,File1.NT(Indices1).Base);
   File2Sequence=cat(2,File2.NT(Indices2).Base);
   
   if abs(length(FASTA(1).Sequence)-length(File1Sequence))>0 || abs(length(FASTA(2).Sequence)-length(File2Sequence))>0
     ErrorMsg='ERROR: Seed alignment does not correspond with pdb files.';
     save(fullfile(ResultFolder, [Query.Name '.mat']), 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg','Query');
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
            align1(end+1)=ct1;
            align2(end+1)=ct2;
         end
      end
   end
elseif ~exist('seed1') || isempty(seed1) %#ok<EXIST>
   [align1 align2] = drNeedlemanWunsch(File1,Indices1,File2,Indices2,.999,2);
else
   for i=1:length(seed1)
      if iscell(seed1)
         tmp=zIndexLookup(File1,seed1{i,1},seed1{i,2},0);
         align1(i)=find(Indices1==tmp);
         tmp=zIndexLookup(File2,seed2{i,1},seed2{i,2},0);
         align2(i)=find(Indices2==tmp);
      elseif isvector(seed1)
         align1=seed1;
         align2=seed2;
      end
   end
end

maxdist=15;
getMoreNeigh = true;
while getMoreNeigh == true;  
  c = cat(1,File1.NT(Indices1).Center);           % nucleotide centers
  File1.Distance = full(zMutualDistance(c,Inf));  % Distance matrix for A
  d = cat(1,File2.NT(Indices2).Center);           % nucleotide centers
  File2.Distance = full(zMutualDistance(d,Inf));  % Distance matrix for B
  A = triu(File1.Distance);                       % Distance matrix is symmetrical;
  B = triu(File2.Distance);                       % don't want duplicate entries returned in next lines
  [rA cA] = find(A<maxdist & A>0);
  [rB cB] = find(B<maxdist & B>0);
  [S1 IX] = sort(rA);
  rA = S1;           
  cA = cA(IX);                                    
  length(rA);
  [S2 IX] = sort(rB);
  rB = S2;
  cB = cB(IX);
  ATrips = DoublesToTriples([rA cA]);
  BTrips = DoublesToTriples([rB cB]);
  AQuads = TriplesToQuads(ATrips);
  BQuads = TriplesToQuads(BTrips);
  clear rA; clear cA; clear rB; clear cB; clear ATrips; clear BTrips;
  if isempty(AQuads) || isempty(BQuads)
     maxdist=maxdist+10;
     getMoreNeigh = true;
     clear AQuads; clear BQuads;
  else
     getMoreNeigh = false;
  end
end

%Balance number of neighborhoods containing each nucleotide
if numNeigh > 0
   AQuads = rGetQuads(AQuads,A,numNeigh);
   BQuads = rGetQuads(BQuads,B,numNeigh);
end
% For each nt in A aligned with a gap, the nearest aligned nt in B is found
% to use as the center of the band for that nt
for i =1:length(Indices1)
   if ~any(align1==i)
      align1 = [align1 i];
      align2 = [align2 0];
   end
end
[S1 IX] = sort(align1);
align1 = S1;
align2 = align2(IX);
  
GappedAs=find(align2==0);
for p=GappedAs
   tmp = 0;
   ct = 0;
   while tmp == 0
      if tmp == 0
         tmp = align2(max(p-ct,1));
      end
      if tmp == 0
         tmp = align2(min(p+ct,length(align2)));
      end
      ct=ct+1;
   end
   align2(p)=tmp;
end

SB=[];
VMI=[];         %This will be a matrix containing the vertices that are created
                  %eg. if VMI(1) =[2 3], there was a vertex created for the 
                  %2nd A neighborhood and 3rd B neighborhood
nA=length(AQuads(:,1));
nB=length(BQuads(:,1));

List=cell(length(A),length(B));     %List{i,j} is a vector indicating the vertices that align
                                    %nucleotide i of A with nucleotide j of B
                                    
% Cutoff values for screening criterion
cut2=4*4*2*discCut*discCut;
cut3=4*4*3*discCut*discCut;
cut4=4*4*4*discCut*discCut;
SA=round(bandwidth/2);
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
            if disc < discCut
               VMI=vertcat(VMI,[i j]);
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
   AlignedNTs1=[];
   AlignedNTs2=[];
   AlignedIndices1=[];
   AlignedIndices2=[];
   ErrorMsg='None aligned';
   save([webroot '/Results' filesep Query.Name filesep Query.Name '.mat'], 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg', 'Query');
   rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
   rWriteAlignmentFasta(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
   return;
end

[EM]=rMakeEdgeMatrix(VMI,List);
clear List;
clear A;
clear B;

maximalClique = rGetMaximalClique(EM);
CLB = length(maximalClique);
    
refineMore=true;
numLoop=0;
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
   VMI(D,:)=[];
   EM(D,:)=[];
   EM(:,D)=[];
end

if strcmpi(cliqueMethod,'full')
    MC = rCliqueBB(EM);
elseif strcmpi(cliqueMethod,'greedy')
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
   AlignedNTs2{i,1}=File2.NT(AlignedIndices2(i)).Number;
   AlignedNTs2{i,2}=File2.NT(AlignedIndices2(i)).Chain;
   AlignedNTs2{i,3}=File2.NT(AlignedIndices2(i)).Base;
end
if length(AlignedIndices1)>4
  rBarDiagram(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,Query.Name);
end
rWriteAlignmentMatrix(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
rWriteAlignmentFasta(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,NTList,Query.Name);
VP.Write=1;
rSuperimposeNucleotides(File1,[AlignedIndices1 setdiff(Indices1,AlignedIndices1)],File2,[AlignedIndices2 setdiff(Indices2,AlignedIndices2)],VP,length(AlignedIndices1),Query.Name);
rAlignmentSpreadsheet(File1,Indices1,File2,Indices2,AlignedIndices1,AlignedIndices2,Query.Name);
% movefile([fr3dpath filesep Query.Name '.pdf'], [webroot filesep 'Results' filesep Query.Name]);
% movefile([fr3dpath filesep Query.Name '.csv'], [webroot filesep 'Results' filesep Query.Name]);
% movefile([fr3dpath filesep Query.Name '.jpg'], [webroot filesep 'Results' filesep Query.Name]);
% movefile([fr3dpath filesep Query.Name '.fasta'], [webroot filesep 'Results' filesep Query.Name]);
% movefile([fr3dpath filesep Query.Name '.txt'], [webroot filesep 'Results' filesep Query.Name]);
% movefile([fr3dpath filesep Query.Name '.pdb'], [webroot filesep 'Results' filesep Query.Name]);
% 
% save(fullfile(ResultFolder, [Query.Name '.mat']), 'AlignedNTs1', 'AlignedNTs2', 'ErrorMsg', 'Query');

% catch ME
% ME.message
% 'error';
% return;
% end