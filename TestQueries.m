% Test for case when no nucleotides are aligned.
Query.Type = 'local';
pdb1 = '4IG8';
pdb2 = '4FY3';
Chain1{1} = 'B';
Nts1{1}   = 'all';
Chain2{1} = '2';
Nts2{1}   = 'all';
Query.SeedName = '';
Disc{1}     = 0.4;
NeighMin{1} = 1;
Band{1}     = 200;
CliqMeth{1} = 'greedy';
Disc{2}     = 0.5;
NeighMin{2} = 3;
Band{2}     = 70;
CliqMeth{2} = 'greedy';
Disc{3}     = 0.5;
NeighMin{3} = 9;
Band{3}     = 20;
CliqMeth{3} = 'greedy';
R3DAlign(pdb1,Chain1,Nts1, pdb2,Chain2,Nts2, Disc,NeighMin,Band,CliqMeth,Query);