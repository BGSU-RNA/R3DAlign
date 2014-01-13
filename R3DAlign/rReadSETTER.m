% rReadSETTER reads the indicated SETTER Al output text file and stores the records in an array.  Each element has three fields:
% Data(n).Header is the header for the nth sequence
% Data(n).Aligned is the sequence, with gaps
% Data(n).Sequence is the sequence with all gaps stripped out

% FASTA =rReadSETTER('ccac5e00-674a-11e3-812a-005056020102.txt');
% FASTA = zReadFASTA('Als\Greengenes_Ecoli_Tth_16S_2_sequences.fasta');

function [AlignedNTs1 AlignedNTs2] = rReadSETTER(Filename)

fid = fopen(Filename,'r');
tline=[];
while isempty(strfind(tline, 'ID'))
    tline = fgetl(fid);
    if ~ischar(tline)
       disp('Problem -1!');
       return; 
    end
end
tline = strrep(tline, 'ID', '');
Filename1 = strrep(tline, ' ', '');
while isempty(strfind(tline, 'chain'))
    tline = fgetl(fid);
    if ~ischar(tline)
       disp('Problem -1!');
       return; 
    end
end
tline = strrep(tline, 'chain', '');
Chain1 = strrep(tline, ' ', '');
while isempty(strfind(tline, 'ID'))
    tline = fgetl(fid);
    if ~ischar(tline)
       disp('Problem -1!');
       return; 
    end
end
tline = strrep(tline, 'ID', '');
Filename2 = strrep(tline, ' ', '');
while isempty(strfind(tline, 'chain'))
    tline = fgetl(fid);
    if ~ischar(tline)
       disp('Problem -1!');
       return; 
    end
end
tline = strrep(tline, 'chain', '');
Chain2 = strrep(tline, ' ', '');

File1 = zAddNTData(Filename1,0);
File2 = zAddNTData(Filename2,0);

c=cat(2,File1.NT.Chain);
Indices1 = find(lower(c)==lower(Chain1));
c=cat(2,File2.NT.Chain);
Indices2 = find(lower(c)==lower(Chain2));

while isempty(strfind(tline, 'Aligned residues:'))
    tline = fgetl(fid);
    if ~ischar(tline)
       disp('Problem 0!');
       return; 
    end
end
while isempty(strfind(tline, 'Aligned residues:'))
    tline = fgetl(fid);
    if ~ischar(tline)
       disp('Problem 1!');
       return; 
    end
end
AlignedNTs1=cell(1);
AlignedNTs2=cell(1);
ct=0;
while 1
    ct=ct+1;
    L = fgetl(fid);
    if L == -1
        break;
    end
    if ~ischar(L)
       disp('Problem 2!');
       return;
    end
    L = strrep(L, '^', '');   %SETTER output is 190^L instead of 190L for NT #'s
    L = strtrim(L);
    k = findstr(L, '-');
    if length(k) ~= 5
       disp('Problem 3!');
       return;
    end
    AlignedNTs1{ct,1} = L(k(1)+1:k(2)-2);     %1st NT #
    AlignedNTs1{ct,2} = L(1:k(1)-1);          %1st Chain
    AlignedNTs1{ct,4} = zIndexLookup(File1,AlignedNTs1{ct,1},AlignedNTs1{ct,2}); %1st index
    
    AlignedNTs2{ct,1} = L(k(5)+1:end);
    AlignedNTs2{ct,2} = L(k(4)+2:k(5)-1);
    AlignedNTs2{ct,4} = zIndexLookup(File2,AlignedNTs2{ct,1},AlignedNTs2{ct,2}); %2nd index
    
    if isempty(AlignedNTs1{ct,4}) || isempty(AlignedNTs2{ct,4})
       AlignedNTs1(ct,:) = [];
       AlignedNTs2(ct,:) = [];
       ct = ct - 1;
    else
       AlignedNTs1{ct,3} = File1.NT(AlignedNTs1{ct,4}).Base;                      %1st Base 
       AlignedNTs2{ct,3} = File2.NT(AlignedNTs2{ct,4}).Base;                      %2nd Base
    end
end
fclose(fid);
Name=[Filename1 '(' Chain1 ')-' Filename2 '(' Chain2 ')_SETTER'];
rAlignmentSpreadsheet(File1,Indices1,File2,Indices2,cell2mat(AlignedNTs1(:,4)),cell2mat(AlignedNTs2(:,4)),Name,'');
rBarDiagram(File1,Indices1,File2,Indices2,cell2mat(AlignedNTs1(:,4)),cell2mat(AlignedNTs2(:,4)),Name,'SETTER');

%SETTER has crossing alignments - remove them
A1=cell2mat(AlignedNTs1(:,4));
A2=cell2mat(AlignedNTs2(:,4));
A1=A1';
A2=A2';

if ~issorted(A1)      %1st set of indices should be in increasing order
    disp('Problem 10!')
end


while ~issorted(A2)
   [i j] = max(abs(tiedrank(A2)-tiedrank(A1)));  %Think of bar diagram - is computing which line crosses the most # of other lines
   A1(j)=[];
   A2(j)=[];
   AlignedNTs1(j,:)=[];
   AlignedNTs2(j,:)=[];
end

ShortOutFilename=[Filename1 '(' Chain1 ')all--' Filename2 '(' Chain2 ')all'];
OutFilename=[Filename1 '(' Chain1 ')all--' Filename2 '(' Chain2 ')all_SETTER_Modified'];
rAlignmentSpreadsheet(File1,Indices1,File2,Indices2,cell2mat(AlignedNTs1(:,4)),cell2mat(AlignedNTs2(:,4)),OutFilename,'');
rBarDiagram(File1,Indices1,File2,Indices2,cell2mat(AlignedNTs1(:,4)),cell2mat(AlignedNTs2(:,4)),OutFilename,'SETTER');
movefile([pwd filesep OutFilename '*'], fullfile(pwd, 'R3D Align Output', OutFilename));
Query.Time=0;
rAnalyzeAlignmentNew(File1,File2,Indices1,Indices2,cell2mat(AlignedNTs1(:,4)),cell2mat(AlignedNTs2(:,4)),OutFilename,ShortOutFilename,'',Query);

return;