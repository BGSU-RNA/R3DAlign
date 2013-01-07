% zCompareAlignmentVectors(Indices1,a1,Indices2,a2) provides diagnostics for alignments given in vectors a1 and a2.  
% These are cell arrays of vectors which tell which indices align.  
% Because indices are sometimes out of linear order in a PDB file, we need the mapping from position in the RNA to the indices, 
% which are given by Indices1 and Indices2.

%zCompareAlignment can be run alternatively by passing a cell array of up to 4 saved .mat alignment filenames
%For example, zrCompareAlignmentVectors({'3V2F(A)all--1S72(0)all_d03_p1_B500_g_d04_p3_B120_g_d05_p9_B20_g','3V2F(A)all--1S72(0)all_d03_p1_B500_g.mat'});
function [average max1] = zCompareAlignmentVectors(Indices1,a1,Indices2,a2,names,m)

if nargin < 2
   Filenames=Indices1;
   for i=1:length(Filenames)
      load(Filenames{i});
%       a1{i} = cat(2,AlignedNTs1{:,4}); %#ok<*USENS>
%       a2{i} = cat(2,AlignedNTs2{:,4});
      a1{i} = AlignedIndices1; %#ok<*USENS>
      a2{i} = AlignedIndices2;
      names{i}=Filenames{i};
      
   end
end

if nargin < 6,
  m = 200;
end

IndexToPos1(Indices1) = 1:length(Indices1);
IndexToPos2(Indices2) = 1:length(Indices2);

figure(6)
clf
figure(7)
clf

T = '';
co = {'red','blue','green','black','magenta','cyan'};

for r = 1:length(a1),
  b1 = IndexToPos1(a1{r});
  b2 = IndexToPos2(a2{r});

  [c1,c2] = rAugmentAlignment(Indices1,b1,b2);

  [y,i] = sort(c1);
  d1{r} = c1(i);
  d2{r} = c2(i); %#ok<*AGROW>

[length(c1) length(c2) length(i) length(d1{r}) length(d2{r})]

  T = [T names{r} ' is ' co{r} ', ']; %#ok<AGROW>

  colors = 'rbgkmc';

  figure(6)
  plot(d1{r},d2{r},colors(min(6,r)));
  hold on

  if r > 1 && r < 6,
    figure(7)
    subplot(2,2,r-1)
    a = d2{r} - d2{1};
    a = min(a,m);
    a = max(a,-m);
    n = hist(a,30);
    hist(a,30);
    axis([-m m 0 1.1*max(n)]);
    if nargin > 4,
      title([names{r} ' versus ' names{1}]);
    end
    average=mean(abs(a));
    max1=max(abs(a));
  end
end
figure(6)
xlabel(T)
