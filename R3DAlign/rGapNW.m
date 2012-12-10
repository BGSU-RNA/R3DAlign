% rGapNW(File1,Indices1,File2,Indices2,p,d,e) aligns sequences using probability p 
% of base conservation and gap open penalty d and gap extension penalty e;

% align1 and align2 list which indices are aligned 
% charAlign1 and charAlign2 give the character alignment, eg.  ACA--GA (charAlign1)
%                                                              AC-UUGA (charAlign2)

function [align1,align2,charAlign1,charAlign2] = rGapNW(File1,Indices1,File2,Indices2,p,d,e)

seq1 = cat(2,File1.NT(Indices1).Base);
seq2 = cat(2,File2.NT(Indices2).Base);

N = length(seq1);
M = length(seq2);

align1 = [];
align2 = [];

pa=1/N;
pb=1/M;

nwmatrix = zeros(N+1,M+1);
tracematrix = zeros(N+1,M+1,'uint8');
GEMatrix = zeros(N+1,M+1,'uint8');      %Gap extension matrix (= 1 if gap has been opened)
GEMatrix(1,2:M+1)=1;
GEMatrix(2:N+1)=1;

for i = 1:N,
  nwmatrix(i+1,1) = -d - i*(e-2);
  tracematrix(i+1,1) = 2;
end

for i = 1:M,
  nwmatrix(1,i+1) = -d - i*(e-2);
  tracematrix(1,i+1) = 1;
end

S1 = seq1'*ones(1,M);
S2 = ones(N,1)*seq2;
S = log(p)*(S1==S2) + log((1-p))*(S1~=S2) - log(pa*pb); %Scoring Matrix

for i = 1:N,
  for j = 1:M,
    x1 = nwmatrix(i+1,j)-d*(GEMatrix(i+1,j)==0)-e*(GEMatrix(i+1,j)==1);
    x2 = nwmatrix(i,j+1)-d*(GEMatrix(i,j+1)==0)-e*(GEMatrix(i,j+1)==1);
    x3 = nwmatrix(i,j)+S(i,j);
%     [a,b] = max([nwmatrix(i+1,j)-d nwmatrix(i,j+1)-d nwmatrix(i,j)+S(i,j)]);
    [a,b] = max([x1 x2 x3]);
    nwmatrix(i+1,j+1) = a;                % max value
    tracematrix(i+1,j+1) = b;             % which direction the max is from
    if b==1 || b==2
       GEMatrix(i+1,j+1)=1;
    end
  end
end

i = N;     % current location in seq1
j = M;     % current location in seq2
matches = 0;
ct = 0;
charAlign1=[];
charAlign2=[];
while (i > 0) || (j > 0),
  ct=ct+1;
  switch tracematrix(i+1,j+1),
  case 1,
    charAlign1=['-' charAlign1]; %#ok<*AGROW>
    charAlign2=[seq2(j) charAlign2];
    j = j - 1;
  case 2,
    charAlign2=['-' charAlign2];
    charAlign1=[seq1(i) charAlign1];
    i = i - 1;
  case 3,
    charAlign1=[seq1(i) charAlign1];
    charAlign2=[seq2(j) charAlign2];
    align1 = [i align1]; %#ok<AGROW>
    align2 = [j align2]; %#ok<AGROW>
    i = i - 1;
    j = j - 1;
    matches = matches + 1;
  end
end

