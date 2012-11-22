% zGetNTData(Filenames,ReadCode,Verbose) reads data files, depending on ReadCode
%
% If ReadCode = 0, it looks for Filename.mat and reads it if it exists.
% If ReadCode = 1, it looks for Filename.mat, reads it, and re-does 
%    the classification of interacting pairs
% If ReadCode = 2, it looks for Filename.mat, reads it 
% If ReadCode = 3, it looks for Filename.mat, reads it, re-does the
%    classification of pairs
% If ReadCode = 4, it reads Filename.pdb, analyzes each nucleotide, 
%    and classifies interacting pairs

function [Files] = zGetNTData(Filenames,ReadCode,Verbose)

[CL, CurrentVersion] = zClassLimits;  % read current version number

if nargin < 2,
  ReadCode = 0;
end

if nargin < 3,
  Verbose = 0;
end

% path(path,pwd);

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};
end

for f=1:length(Filenames),
  Filename = Filenames{f};

  if isempty(strfind(lower(Filename),'.pdb')),
    PDBFilename  = [Filename '.pdb'];
  else
    PDBFilename  = Filename;
    i = strfind(Filename,'.');
    Filename = Filename(1:(i(end)-1));
    ReadCode = 4;
  end

  FILENAME = upper(Filename);
  filename = lower(Filename);

  ClassifyCode = 0;
  SaveCode     = 0;
  ReadPDB      = 0;

  if (ReadCode > 3),
    ReadFull = 1;
  else
    ReadFull = 0;
  end

  clear File

  if ReadCode == 4,                     % re-read the PDB file
    File = zReadandAnalyze(PDBFilename,Verbose);   % might not work on a Mac
    ClassifyCode = 1;
    ReadFull = 0;                       % no need to read PDB file again
    ReadPDB  = 1;
  else                                  % try to load a precomputed version
    if (exist(strcat(Filename,'.mat'),'file') > 0),
      load(strcat(Filename,'.mat'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded %s\n', [Filename '.mat']);
      end
    elseif (exist(strcat(FILENAME,'.MAT'),'file') > 0),  % helps on a Mac
      load(strcat(FILENAME,'.MAT'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded  %s\n', [FILENAME '.MAT']);
      end
    elseif (exist(strcat(filename,'.mat'),'file') > 0),  % helps on a Mac
      load(strcat(filename,'.mat'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded   %s\n', [filename '.MAT']);
      end
    else
      ReadFull = 1;
    end
  end

  if (ReadFull == 1),                     % try to load the full version
    File = zReadandAnalyze(PDBFilename,Verbose);
    ClassifyCode = 1;
    ReadPDB      = 1;
  end

  if ~exist('File'),
    File = [];
  end

  if ~isfield(File,'ClassVersion'),
    File.ClassVersion = 0;
  end

  if ~isfield(File,'BasePhosphate'),
    File.BasePhosphate = sparse(zeros(length(File.NT)));
  end

  if isempty(File.BasePhosphate),
    File.BasePhosphate = sparse(zeros(length(File.NT)));
  end

  if ~isfield(File,'Covalent'),
    if length(File.NT) > 1,
      File = zBackboneContinuity(File);
    else
      File.Covalent = sparse(zeros(length(File.NT)));
    end
  end

  if isfield(File,'Inter'),
    File = rmfield(File,'Inter');
  end

  if isfield(File,'SizeCode'),
    File = rmfield(File,'SizeCode');
  end

  if isfield(File,'Pair'),
    File = rmfield(File,'Pair');
  end

  if isfield(File,'Header'),
    File = rmfield(File,'Header');
  end

  Overlap = 0;

  File.Distance = [];                        % only compute this if needed

  if length(File.NT) > 1,                    % if it has nucleotides,
    if length(File.NT(1).Sugar(:,1)) < 13,
      File = zStoreO3(File);
    end

    if (ReadCode == 1) | (ReadCode == 3) | (ReadCode == 4) | ... 
      (ClassifyCode == 1) | (File.ClassVersion < CurrentVersion),

      File.Edge = sparse(File.NumNT,File.NumNT);
      File.Coplanar = sparse(File.NumNT,File.NumNT);

      c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers

      File.Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
      d = sort(nonzeros(File.Distance));

      if isempty(d),
       % more than one nucleotide, but too far apart to pair
       File.Range = sparse(File.NumNT,File.NumNT);
       File.BasePhosphate = sparse(File.NumNT,File.NumNT);
      else
       if d(min(10,length(d))) < 1,
         fprintf('%s has overlapping nucleotides and should be avoided\n',File.Filename);

         Overlap = 1;
       else
         t = cputime;
         File = zClassifyPairs(File,Verbose);

         if Verbose > 0,
           fprintf(' Base-phosphate interactions ...');
         end

         File = zPhosphateInteractions(File);
         File.ClassVersion = CurrentVersion;
         ClassifyCode = 1;

         if Verbose > 0,
           fprintf(' Interaction range ... ');
         end

         File = zInteractionRange(File,Verbose);

         if Verbose > 0,
           fprintf('\nClassification took %4.2f minutes\n', (cputime-t)/60);
         end
       end
      end
    end
  else
    File.ClassVersion = CurrentVersion;
  end

  if length(File.NT) > 0,
   if ~isfield(File.NT(1),'Syn'),
    SynList = mSynList(File);
    for k=1:length(File.NT),
      File.NT(k).Syn = SynList(k);
    end
    SaveCode = 1;
   end
  end

  if ~isfield(File,'PDBFilename'),
    File.PDBFilename = [File.Filename '.pdb'];   % will self-correct the
                                                 % next time PDB is read
  end

  if ~isfield(File,'Backbone'),
    File = zBackboneConformation(File,Verbose);
    if File.NumNT > 1,
      File.Backbone(1,1) = 1;                    % register that we checked
    end
    SaveCode = 1;
  end

  if File.NumNT > 1 && sum(sum(File.Backbone)) == 0,
    File = zBackboneConformation(File,Verbose);
    SaveCode = 1;
  end

  if ~isfield(File,'Range') || ~isfield(File,'Crossing'),
    File = zInteractionRange(File,Verbose);
    ClassifyCode = 1;
  end

  if ~isfield(File,'Flank'),
    File = xFlankingPairs(File);
  end

  if ~isfield(File,'Info'),
    File = zGetPDBInfo(File);          % get resolution and other info
    SaveCode = 1;
  else
    if isempty(File.Info.Descriptor),
      File = zGetPDBInfo(File);          % look for file information
      if ~isempty(File.Info.Descriptor),
        SaveCode = 1;
      end
    end
  end

  % ----------------- If it just read the .pdb file and no pairs, look at BUC

  if (ReadPDB == 1) && (isempty(strfind(PDBFilename,'.pdb1'))),
    E = abs(File.Edge);                       % all interactions
    bp = full(sum(sum((E > 0) .* (E < 15))));       % number of basepairs
    r  = bp / length(File.NT);                % ratio of bp to nt
    if (r < 0.4),                             % very few basepairs found
      if Verbose > 0,
        fprintf('Few basepairs found (ratio %7.2f), reading the biological unit coordinates\n',r);
      end
      File1 = zGetNTData([PDBFilename '1'],4,1);   % read biological unit coords
      if length(File1.NT) > 0,
        E1  = abs(File1.Edge);
        bp1 = full(sum(sum((E > 0) .* (E < 15))));       % number of basepairs
        r1  = bp / length(File1.NT);               % ratio of bp to nt
        if Verbose > 0,
          fprintf('Biological unit coordinates have ratio %7.2f\n',r1);
        end
        if r1 > r,
          File = File1;                         % use biological unit coords
        end
      end
    end
  end
    
  File = orderfields(File);

  Saved = 0;

  if (ReadCode > 0) || (ClassifyCode > 0) || (SaveCode > 0),
    zSaveNTData(File,Verbose);
    Saved = 1;
  end

  Files(f) = File;

end
