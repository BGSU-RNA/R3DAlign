function [BPcount,BPnearcount] = rInteractionSummary(File,NTList)
if ischar(File),
  fprintf('Loading PDB Info...\n');
  Filename = File;
  File = zGetNTData(Filename,0);
end

% if NTList is a cell array of numbers, look up the indices
if ischar(NTList)
   if isequal('all',lower(strtrim(NTList)))
      NTList=1:File.NumNT;
   elseif length(NTList)==1
      chain=NTList;
      c=cat(2,File.NT.Chain);
      NTList = find(lower(c)==lower(chain));
      if isempty(NTList)
         fprintf('Invalid chain entered for the first structure.\n');
         return;
      end
   else
      NTList = {NTList};
   end  
end

if iscell(NTList),
   [Indices chains] = zIndexLookup(File,NTList);
   for i=1:length(chains)
      if length(chains{i})>1
         Chain1 = input('Enter the chain identifier for the first structure: ','s');
         Indices = zIndexLookup(File,NTList,Chain1);
         break;
      end
   end
else
  Indices = NTList;
end  
E  = fix(abs(File.Edge));
E=E(Indices,Indices);
E=triu(E);
B  = E .* (E > 0) .* (E < 24);   % pairs and stacks
C  = File.Crossing;
C=C(Indices,Indices);

Names={'cww','tww','cwh','twh','cws','tws','chh','thh','chs','ths','css','tss'};
BPcount=cell(12,3);
BPnearcount=cell(12,3);
for i=1:12
    BPcount{i,1}=length(find(E==i));                        %True Pairs
%     BPcount{i,2}=length(find((E==i).*(C==0)));               %nested pairs
%     BPcount{i,3}=round2(BPcount{i,2}/BPcount{i,1},.01);
%     BPcount{i,3}=length(find((E==i).*(C>0)));              %nonnested pairs
    BPcount{i,2}=length(find((E==i).*(C>3)));               %long range pairs
    BPcount{i,3}=round2(BPcount{i,2}/BPcount{i,1},.01);
    BPnearcount{i,1}=length(find(E==(100+i)))+BPcount{i,1};              %near Pairs
%     BPnearcount{i,2}=length(find((E==(100+i)).*(C==0)))+BPcount{i,2};     %nested near pairs
%     BPnearcount{i,3}=round2(BPnearcount{i,2}/BPnearcount{i,1},.01);
%     BPnearcount{i,3}=length(find((E==(100+i)).*(C>0)))+BPcount{i,3};    %nonnested near pairs
    BPnearcount{i,2}=length(find((E==(100+i)).*(C>3)))+BPcount{i,2};           %long range near pairs
    BPnearcount{i,3}=round2(BPnearcount{i,2}/BPnearcount{i,1},.01);
    
end
for i=1:12
   BPFamily(i).Name=Names{i}; %#ok<*AGROW,*AGROW>
   BPFamily(i).freqs=zeros(5,5);
   BPFamily(i).relfreqs=zeros(5,5); %#ok<*AGROW>
   BPFamily(i).nearfreqs=zeros(5,5);
   BPFamily(i).nearrelfreqs=zeros(5,5);
   BPFamily(i).freqsLR=zeros(5,5);
   BPFamily(i).nearfreqsLR=zeros(5,5);
end
P  = E .* (E > 0) .* (E < 13);   % pairs
[i j k] = find(P);
for m=1:length(i)
   a1=File.NT(Indices(i(m))).Base;
   b1=File.NT(Indices(j(m))).Base;
   a1=char2num(a1);
   b1=char2num(b1);
   BPFamily(k(m)).freqs(a1,b1)=BPFamily(k(m)).freqs(a1,b1)+1;
   BPFamily(k(m)).freqs(a1,5)=BPFamily(k(m)).freqs(a1,5)+1;
   BPFamily(k(m)).freqs(5,b1)=BPFamily(k(m)).freqs(5,b1)+1;
   BPFamily(k(m)).freqs(5,5)=BPFamily(k(m)).freqs(5,5)+1;
end

for i=1:12
   BPFamily(i).nearfreqs=BPFamily(i).freqs;
end
NP  = E .* (E > 100) .* (E < 113);   % pairs
[i j k] = find(NP);
for m=1:length(i)
   a1=File.NT(Indices(i(m))).Base;
   b1=File.NT(Indices(j(m))).Base;
   a1=char2num(a1);
   b1=char2num(b1);
   BPFamily(k(m)-100).nearfreqs(a1,b1)=BPFamily(k(m)-100).nearfreqs(a1,b1)+1;
   BPFamily(k(m)-100).nearfreqs(a1,5)=BPFamily(k(m)-100).nearfreqs(a1,5)+1;
   BPFamily(k(m)-100).nearfreqs(5,b1)=BPFamily(k(m)-100).nearfreqs(5,b1)+1;
   BPFamily(k(m)-100).nearfreqs(5,5)=BPFamily(k(m)-100).nearfreqs(5,5)+1;
end
Label{1,2}='A';Label{1,3}='C';Label{1,4}='G';Label{1,5}='U';
Label{2,1}='A';Label{3,1}='C';Label{4,1}='G';Label{5,1}='U';
Label{5,5}='';
for i=1:12
   Label{1,1}=Names{i};
   if mod(i,2)==1
      xlswrite([File.Filename '_Interaction Summary'],Label,'Freqs',['A' num2str((i-1)/2*7+1) ':E' num2str((i-1)/2*7+5)]);
      xlswrite([File.Filename '_Interaction Summary'],BPFamily(i).freqs,'Freqs',['B' num2str((i-1)/2*7+2) ':F' num2str((i-1)/2*7+6)]);
      xlswrite([File.Filename '_Interaction Summary'],Label,'NearFreqs',['A' num2str((i-1)/2*7+1) ':E' num2str((i-1)/2*7+5)]);
      xlswrite([File.Filename '_Interaction Summary'],BPFamily(i).nearfreqs,'NearFreqs',['B' num2str((i-1)/2*7+2) ':F' num2str((i-1)/2*7+6)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'RelFreqs',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).relfreqs,'RelFreqs',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'RelFreqsNear',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).nearrelfreqs,'RelFreqsNear',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'FreqsLR',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).freqsLR,'FreqsLR',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'FreqsLRNear',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).nearfreqsLR,'FreqsLRNear',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
   else
      xlswrite([File.Filename '_Interaction Summary'],Label,'Freqs',['H' num2str((i-2)/2*7+1) ':L' num2str((i-2)/2*7+5)]);
      xlswrite([File.Filename '_Interaction Summary'],BPFamily(i).freqs,'Freqs',['I' num2str((i-2)/2*7+2) ':M' num2str((i-2)/2*7+6)]);
      xlswrite([File.Filename '_Interaction Summary'],Label,'NearFreqs',['H' num2str((i-2)/2*7+1) ':L' num2str((i-2)/2*7+5)]);
      xlswrite([File.Filename '_Interaction Summary'],BPFamily(i).nearfreqs,'NearFreqs',['I' num2str((i-2)/2*7+2) ':M' num2str((i-2)/2*7+6)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'RelFreqs',['T' num2str((i-2)/2*19+1) ':AJ' num2str((i-2)/2*19+17)]);
%       xlswrite([File1.Filename '_' File2.Filename
%       '_alignment_data.xlsx'],BPFamily(i).relfreqs,'RelFreqs',['U' num2str((i-2)/2*19+2) ':AK' num2str((i-2)/2*19+18)]); 
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'RelFreqsNear',['T' num2str((i-2)/2*19+1) ':AJ' num2str((i-2)/2*19+17)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).nearrelfreqs,'RelFreqsNear',['U' num2str((i-2)/2*19+2) ':AK' num2str((i-2)/2*19+18)]); 
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'FreqsLR',['T' num2str((i-2)/2*19+1) ':AJ' num2str((i-2)/2*19+17)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).freqsLR,'FreqsLR',['U' num2str((i-2)/2*19+2) ':AK' num2str((i-2)/2*19+18)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'FreqsLRNear',['T' num2str((i-2)/2*19+1) ':AJ' num2str((i-2)/2*19+17)]);
%       xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).nearfreqsLR,'FreqsLRNear',['U' num2str((i-2)/2*19+2) ':AK' num2str((i-2)/2*19+18)]);
   end
end

Color = (B==1).*(C==0) + 2*(B>1).*(B<13).*(C==0) + 3*(B==1).*(C>0) + 4*(B > 1).*(B < 13) .*(C>0) + 5*(B > 20) .* (B < 25);
                                        % disjoint possibilities

BP = abs(File.BasePhosphate);           % 
BP = (BP > 0) .* (BP < 100);            % exclude near BP and self interactions

[i,j,c] = find(triu(Color));

[ii,jj,cc] = find(6*BP);                % handle this separately, don't add

%[ii,jj,cc] = find(6*BP.*(R > 1));     % range 0 BPh only

k = find(ii ~= jj);                     % eliminate self interactions

% i = [i; ii(k)];
% j = [j; jj(k)];
c = [c; cc(k)];

% cww = length(find(c == 1))
% Tally(1,1) = cww;
% noncww = length(find(c == 2))
% Tally(1,2) = noncww;
% nonne�M�K)6�rNUW�*�Y��q����OI��w�����1��u6��S{�i֭Ye<�f�\[��8l�+{vy�5/�|��p�D�YQ_��1�X���V�2�I��k����y�g��ٌ$���3�xqh�9�d*���ڴ�s%.22H/�������7��*1�t� ��67��Bm�$NŁ�����l$3�
!ᕙr��B3�x�Q�9J�X/���,uT���S9Kt�=����ѷ���://�>p^7�#����m���N7v�v�~̌�)шZ+-4e)|Y�,`��
d�����[#B�i�{��J�8B.��o_*/��-ػ�e��R;��2?�=if�������.�lrq��$P� 3U�#�����Ƨs����Da氚@cY��	H�\��%Lwn��(c%�*�W;�/�t�^<7�e�&@��C�3�2bQڵ�@�*M���̹�&#��z���^Ta+�G%����I{�+,_L�.�.�����r���y��]�M�K)6�rNUW�*�Y��q����OI��w�����1��u6��S{�i֭Ye<�f�\[��8l�+{vy�5/�|��p�D�YQ_��1�X���V�2�I��k����y�g��ٌ$���3�xqh�9�d*���ڴ�s%.22H/�������7��*1�t� ��67��Bm�$NŁ�����l$3�
!ᕙr��B3�x�Q�9J�X/���,uT���S9Kt�=����ѷ���://�>p^7�#����m���N7v�v�~̌�)шZ+-4e)|Y�,`��
d�����[#B�i�{��J�8B.��o_*/��-ػ�e��R;��2?�=if�������.�lrq��$P� 3U�#�����Ƨs����Da氚@cY��	H�\��%Lwn��(c%�*�W;�/�t�^<7�e�&@��C�3�2bQڵ�@�*M���̹�&#��z���^Ta+�G%����I{�+,_L�.�.�����r���y��]֓ݑ�=4�1�y��8��&������n��M�=��q�����H�}����C�0���� ��	�a�%:��~1'̧���΂�D6��4�RXבg?ӭ0ȯ�|�	
A�Ұ�ė�d�]%�%�`�ncx-���.Z��f�T�����N��6= f0/���j�L-�t -|�K�d�66��'�
���ab��C��n�Ŋ(h�=�>85��<��9�ɹደ����������ګ���i�� O]f�Z4oΗn>���^����������� 8�t4N�y��4��%�"ͥwQ~ �&)#�P�c���w�@Ƈވ�ܙ�K���9��P*�'?