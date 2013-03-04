% Test for case when no nucleotides are aligned.
clear
Query.Name = '5134bd88ed961';
Query.Type = 'web';
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
webWrapper(pdb1,Chain1,Nts1, pdb2,Chain2,Nts2, Disc,NeighMin,Band,CliqMeth,Query);

%Molecule 2 isn't read in correctly
clear
fprintf('\n');
Query.Name = '5134c80367b92';
Query.Type = 'web';
pdb1  = 'uploaded';
Name1 = '5134c80367b92_1.pdb';
Query.UploadName1 = '5134c80367b92_1';
pdb2 = 'uploaded';
Name2 = '5134c80367b92_2.pdb';
Query.UploadName2 = '5134c80367b92_2';
Chain1{1} = 'A';
Nts1{1}   = 'all';
Chain2{1} = 'A';
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
webWrapper(pdb1,Chain1,Nts1, pdb2,Chain2,Nts2, Disc,NeighMin,Band,CliqMeth,Query);

%Molecule 1 read in incorrectly
clear
Query.Name = '5134bdfcdb5b3';
Query.Type = 'web';
pdb1  = 'uploaded';
Name1 = '5134bdfcdb5b3_1.pdb';
Query.UploadName1 = '5134bdfcdb5b3_1';
pdb2 = 'uploaded';
Name2 = '5134bdfcdb5b3_2.pdb';
Query.UploadName2 = '5134bdfcdb5b3_2';
Chain1{1} = 'all';
Nts1{1}   = 'all';
Chain2{1} = 'all';
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
webWrapper(pdb1,Chain1,Nts1, pdb2,Chain2,Nts2, Disc,NeighMin,Band,CliqMeth,Query);