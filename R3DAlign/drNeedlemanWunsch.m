% drNeedlemanWunsch(seq1,seq2,p,d) aligns sequences using probability p of base conservation and gap penalty d.  
% align1 and align2 list which indices are aligned  
function [align1,align2] = drNeedlemanWunsch(File1,Indices1,File2,Indices2,p,d)

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

for i = 1:N,
  nwmatrix(i+1,1) = -i*d;
  tracematrix(i+1,1) = 2;
end

for i = 1:M,
  nwmatrix(1,i+1) = -i*d;
  tracematrix(1,i+1) = 1;
end

S1 = seq1'*ones(1,M);
S2 = ones(N,1)*seq2;
S = log(p)*(S1==S2) + log((1-p))*(S1~=S2) - log(pa*pb); %Scoring Matrix

for i = 1:N,
  for j = 1:M,
    [a,b] = max([nwmatrix(i+1,j)-d nwmatrix(i,j+1)-d nwmatrix(i,j)+S(i,j)]);
    nwmatrix(i+1,j+1) = a;                % max value
    tracematrix(i+1,j+1) = b;             % which direction the max is from
  end
end

i = N;     % current location in seq1
j = M;     % current location in seq2
matches = 0;
ct = 0;

while (i > 0) || (j > 0),
  ct=ct+1;
  switch tracematrix(i+1,j+1),
  case 1,
    j = j - 1;
  case 2,
    i = i - 1;
  case 3,
    align1 = [i align1]; %#ok<AGROW>
    align2 = [j align2]; %#ok<AGROW>
    i = i - 1;
    j = j - 1;
    matches = matches + 1;
  end
end