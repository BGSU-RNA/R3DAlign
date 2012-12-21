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
% nonneİMıK)6ÚrNUW*êœYˆ‰qÆĞ÷ïOI¨éwŞêÀ´Å1ßÜu6˜şS{ÇiÖ­Ye<²fü\[Ñë8lé+{vyâ5/”|Š¬p’DîYQ_Æñ1áXÓ¬VÁ2üIş«kşéçĞyôgÊéÙŒ$œŸ¦3§xqh‚9³d*ª‘ÓÚ´Ús%.22H/µº„»¸éºÆ7¥Ç*1ËtÓ ü•67ÇæBm†$NÅ¤“¾ÛÊl$3¦
!á•™rßèB3îx—QÑ9JƒX/Ùğ,uT¹–¢S9Ktø=®™‹ŸÑ·µ¤­://‹>p^7»î«¦#Œ”çĞmŞõ‚N7vövÃ~ÌŒ€)ÑˆZ+-4e)|Y,`ñÃ
d®²ªú[#BİiÉï£—{›ÌJ£8B.õào_*/â“ã-Ø»úe§è§R;‚2?Ÿ=ifÖ®ÑàœÃ¯.·lrqºØ$Pì 3U#ÿœµìÆ§sı“•úDaæ°š@cYëÃ	H\õƒ%Lwn‘(c%ã*ÂW;ô/¡t…^<7¼eú&@ìCÔ3Ä2bQÚµ±@¯*MÃû¿Ì¹ˆ&#ÄÿzŒ¿Æ^Ta+¨G%¶ë×ÇI{Ø+,_L–.¨.ˆƒıñ´’¸rà™Âyèäµ]İMıK)6ÚrNUW*êœYˆ‰qÆĞ÷ïOI¨éwŞêÀ´Å1ßÜu6˜şS{ÇiÖ­Ye<²fü\[Ñë8lé+{vyâ5/”|Š¬p’DîYQ_Æñ1áXÓ¬VÁ2üIş«kşéçĞyôgÊéÙŒ$œŸ¦3§xqh‚9³d*ª‘ÓÚ´Ús%.22H/µº„»¸éºÆ7¥Ç*1ËtÓ ü•67ÇæBm†$NÅ¤“¾ÛÊl$3¦
!á•™rßèB3îx—QÑ9JƒX/Ùğ,uT¹–¢S9Ktø=®™‹ŸÑ·µ¤­://‹>p^7»î«¦#Œ”çĞmŞõ‚N7vövÃ~ÌŒ€)ÑˆZ+-4e)|Y,`ñÃ
d®²ªú[#BİiÉï£—{›ÌJ£8B.õào_*/â“ã-Ø»úe§è§R;‚2?Ÿ=ifÖ®ÑàœÃ¯.·lrqºØ$Pì 3U#ÿœµìÆ§sı“•úDaæ°š@cYëÃ	H\õƒ%Lwn‘(c%ã*ÂW;ô/¡t…^<7¼eú&@ìCÔ3Ä2bQÚµ±@¯*MÃû¿Ì¹ˆ&#ÄÿzŒ¿Æ^Ta+¨G%¶ë×ÇI{Ø+,_L–.¨.ˆƒıñ´’¸rà™Âyèäµ]Ö“İ‘Ù=4¨1Êy¶8®ğ&’ÔËÔàõnäMâ=±ÛqóØøìùHü}ßüØCê·0×íØ  ¤	aâ²%:À¾~1'Ì§¡§´Î‚ËD6Ï×4ÕRX×‘g?Ó­0È¯ù|Ÿ	
AÛÒ°¤Ä—¼d‘]%ô%€`Ñncx-Í€å¢.ZÙÚf€T²ôµÎN¸â6= f0/ıÑçjñL-óšt -|´KÎdò66¬'ì
ÔÕıab¨¡C¶¬n»ÅŠ(h©= >85ªç<³‡9ÚÉ¹á‹°åÙîôÃ§÷å¨ÉÚ«˜‰×iµ„ O]fè’Z4oÎ—n>¶ÊÉ^³£Õ÷†âÚãó‹» 8Ùt4N™yöÍ4Ê¢%Ö"Í¥wQ~ ÿ&)#ÎP‘c¡ûwÔ@Æ‡ŞˆáÜ™èKŒèï9¦ÁP*”'?