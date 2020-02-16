# About

Comparing 3D structures of homologous RNA molecules yields information about sequence and structural variability. To compare large RNA 3D structures, accurate automatic comparison tools are needed. In this article, we introduce a new algorithm and web server to align large homologous RNA structures nucleotide by nucleotide using local superpositions that accommodate the flexibility of RNA molecules. Local alignments are merged to form a global alignment by employing a maximum clique algorithm on a specially defined graph that we call the ‘local alignment’ graph.

The algorithm is implemented in a program suite and web server called ‘R3D Align’. The R3D Align alignment of homologous 3D structures of 5S, 16S and 23S rRNA was compared to a high-quality hand alignment. A full comparison of the 16S alignment with the other state-of-the-art methods is also provided. The R3D Align program suite includes new diagnostic tools for the structural evaluation of RNA alignments. The R3D Align alignments were compared to those produced by other programs and were found to be the most accurate, in comparison with a high quality hand-crafted alignment and in conjunction with a series of other diagnostics presented. The number of aligned base pairs as well as measures of geometric similarity are used to evaluate the accuracy of the alignments.

# Installation

From the command line.  Skip the first line if FR3D is already installed.

    git clone git@github.com:BGSU-RNA/FR3D.git
    git clone git@github.com:BGSU-RNA/R3DAlign.git

# Usage

The main program is _R3DAlign.m_. After starting Matlab or Octave:

    cd FR3D;
    addpath FR3DSource PrecomputedData PDBFiles SearchSaveFiles;
    cd R3DAlign;
    addpath(genpath(pwd));
    % test alignment of two 5S rRNAs
    [a1,a2] = R3DAlign('2AW4',{'A'},{'all'},'2J01',{'B'},{'all'},0.5,9,50,'greedy');

Alternatively, manually set the Matlab path to include the FR3D folders listed above.  This will allow R3D Align to use binary .mat files for structures already downloaded and annoted by FR3D.  However, any new .mat files will be stored in the R3DAlign folder, not back in the FR3D folder.

To produce pdb files, use the following commands:

    Query.Type = 'local';
    Query.Name = 'output_file'; % will produce output_file.pdb in the current working directory
    [a1,a2] = R3DAlign('2AW4',{'A'},{'all'},'2J01',{'B'},{'all'},0.5,9,50,'greedy',Query);

# Credits

Developed by Ryan Rahrig.
Transferred to Github by Anton Petrov.
Updated by Craig Zirbel

Tested on Matlab R2007b, R2019b, and Octave 3.6.3, at various times.  We try to maintain backward compatibility, but it cannot be assured.