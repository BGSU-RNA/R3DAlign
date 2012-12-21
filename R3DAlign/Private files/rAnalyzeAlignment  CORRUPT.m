% Evaluate the number of near basepairs aligned with  true basepairs, near with near, 
% true with true, etc.

function [BPFamily, T] = rAnalyzeAlignment(File1,File2,SpreadsheetName)
clear a;
clear b;
[a, b]=xlsread(SpreadsheetName); %#ok<ASGLU>
for k = [2 5]
   for i=1:length(b(:,2))
      b(i,2)=lower(strtrim(b(i,2)));
      b(i,5)=lower(strtrim(b(i,5)));
   end
   
   for i=1:length(b(:,k))
      if isequal(b{i,k},'tsh')
         b{i,k}='ths';
      elseif isequal(b{i,k},'csh')
         b{i,k}='chs';
      elseif isequal(b{i,k},'tsw')
         b{i,k}='tws';
      elseif isequal(b{i,k},'csw')
         b{i,k}='cws';
      elseif isequal(b{i,k},'thw')
         b{i,k}='twh';
      elseif isequal(b{i,k},'chw')
         b{i,k}='cwh';
      elseif isequal(b{i,k},'ntsh')
         b{i,k}='nths';
      elseif isequal(b{i,k},'ncsh')
         b{i,k}='nchs';
      elseif isequal(b{i,k},'ntsw')
         b{i,k}='ntws';
      elseif isequal(b{i,k},'ncsw')
         b{i,k}='ncws';
      elseif isequal(b{i,k},'nthw')
         b{i,k}='ntwh';
      elseif isequal(b{i,k},'nchw')
         b{i,k}='ncwh';
      elseif ~isempty(b{i,k}) && isequal(b{i,k}(1),'s')
         b{i,k}='stacking';
      end
   end
end
Names={'cww','tww','cwh','twh','cws','tws','chh','thh','chs','ths','css','tss'};
for i=1:12
   BPFamily(i).Name=Names{i}; %#ok<*AGROW>
   BPFamily(i).numAligned=0;
   BPFamily(i).freqs=zeros(17,17);
   BPFamily(i).relfreqs=zeros(17,17);
   BPFamily(i).nearfreqs=zeros(17,17);
   BPFamily(i).nearrelfreqs=zeros(17,17);
   BPFamily(i).freqsLR=zeros(17,17);
   BPFamily(i).nearfreqsLR=zeros(17,17);
end

for i=2:length(b(:,2))
   a1=b{i,1}(3);
   b1=b{i,3}(3);
   a2=b{i,4}(3);
   b2=b{i,6}(3);
   a1=char2num(a1);
   b1=char2num(b1);
   a2=char2num(a2);
   b2=char2num(b2);
   if a1~=0 && b1~=0 && zIndexLookup(File1,b{i,1}(4:end),b{i,1}(1)) < zIndexLookup(File1,b{i,3}(4:end),b{i,1}(1))
      for j=1:12
         if isequal(strtrim(b{i,2}),Names{j})
            if isequal(b{i,5},Names{j})
               BPFamily(j).numAligned=BPFamily(j).numAligned+1;
               BPFamily(j).freqs((a1-1)*4+b1,(a2-1)*4+b2)=BPFamily(j).freqs((a1-1)*4+b1,(a2-1)*4+b2)+1;
               BPFamily(j).freqs((a1-1)*4+b1,17)=BPFamily(j).freqs((a1-1)*4+b1,17)+1;
               BPFamily(j).freqs(17,(a2-1)*4+b2)=BPFamily(j).freqs(17,(a2-1)*4+b2)+1;
               BPFamily(j).freqs(17,17)=BPFamily(j).freqs(17,17)+1;
            
               if File1.Crossing(zIndexLookup(File1,b{i,1}(4:end),b{i,1}(1)),zIndexLookup(File1,b{i,3}(4:end),b{i,1}(1))) > 3
                  BPFamily(j).freqsLR((a1-1)*4+b1,(a2-1)*4+b2)=BPFamily(j).freqsLR((a1-1)*4+b1,(a2-1)*4+b2)+1;
                  BPFamily(j).freqsLR((a1-1)*4+b1,17)=BPFamily(j).freqsLR((a1-1)*4+b1,17)+1;
                  BPFamily(j).freqsLR(17,(a2-1)*4+b2)=BPFamily(j).freqsLR(17,(a2-1)*4+b2)+1;
                  BPFamily(j).freqsLR(17,17)=BPFamily(j).freqsLR(17,17)+1;
               end
            end
         end
         if isequal(strtrim(b{i,2}),Names{j}) || isequal(strtrim(b{i,2}),['n' Names{j}]) 
            if isequal(b{i,5},Names{j}) || isequal(strtrim(b{i,5}),['n' Names{j}])
               BPFamily(j).nearfreqs((a1-1)*4+b1,(a2-1)*4+b2)=BPFamily(j).nearfreqs((a1-1)*4+b1,(a2-1)*4+b2)+1;
               BPFamily(j).nearfreqs((a1-1)*4+b1,17)=BPFamily(j).nearfreqs((a1-1)*4+b1,17)+1;
               BPFamily(j).nearfreqs(17,(a2-1)*4+b2)=BPFamily(j).nearfreqs(17,(a2-1)*4+b2)+1;
               BPFamily(j).nearfreqs(17,17)=BPFamily(j).nearfreqs(17,17)+1;
               if File1.Crossing(zIndexLookup(File1,b{i,1}(4:end),b{i,1}(1)),zIndexLookup(File1,b{i,3}(4:end),b{i,1}(1))) > 3
                  BPFamily(j).nearfreqsLR((a1-1)*4+b1,(a2-1)*4+b2)=BPFamily(j).nearfreqsLR((a1-1)*4+b1,(a2-1)*4+b2)+1;
                  BPFamily(j).nearfreqsLR((a1-1)*4+b1,17)=BPFamily(j).nearfreqsLR((a1-1)*4+b1,17)+1;
                  BPFamily(j).nearfreqsLR(17,(a2-1)*4+b2)=BPFamily(j).nearfreqsLR(17,(a2-1)*4+b2)+1;
                  BPFamily(j).nearfreqsLR(17,17)=BPFamily(j).nearfreqsLR(17,17)+1;
               end
            end
         end
      end
   end
end
for k=1:12
   for i=1:16
      for j=1:16
       BPFamily(k).relfreqs(i,j)=round2(BPFamily(k).freqs(i,j)/BPFamily(k).freqs(i,17),.01);
       BPFamily(k).nearrelfreqs(i,j)=round2(BPFamily(k).nearfreqs(i,j)/BPFamily(k).nearfreqs(i,17),.01);
      end
   end
   BPFamily(k).relfreqs(:,17)=BPFamily(k).freqs(:,17);
   BPFamily(k).relfreqs(17,:)=BPFamily(k).freqs(17,:);
   BPFamily(k).nearrelfreqs(:,17)=BPFamily(k).nearfreqs(:,17);
   BPFamily(k).nearrelfreqs(17,:)=BPFamily(k).nearfreqs(17,:);
end       

Label{1,2}='AA';Label{1,3}='AC';Label{1,4}='AG';Label{1,5}='AU';
Label{1,6}='CA';Label{1,7}='CC';Label{1,8}='CG';Label{1,9}='CU';
Label{1,10}='GA';Label{1,11}='GC';Label{1,12}='GG';Label{1,13}='GU';
Label{1,14}='UA';Label{1,15}='UC';Label{1,16}='UG';Label{1,17}='UU';
Label{2,1}='AA';Label{3,1}='AC';Label{4,1}='AG';Label{5,1}='AU';
Label{6,1}='CA';Label{7,1}='CC';Label{8,1}='CG';Label{9,1}='CU';
Label{10,1}='GA';Label{11,1}='GC';Label{12,1}='GG';Label{13,1}='GU';
Label{14,1}='UA';Label{15,1}='UC';Label{16,1}='UG';Label{17,1}='UU';
Label{17,17}='';

Temp=cell(40,40);
Temp(3:14,1)=Names';
Temp{1,2}='TRUE';
Temp{2,2}='All';
Temp{2,3}='Long Range';
Temp(3:14,5)=Names';
Temp{1,6}='TRUE and NEAR';
Temp{2,6}='All';
Temp{2,7}='Long Range';
for i=1:12
    Temp{i+2,2}=BPFamily(i).freqs(17,17);
    Temp{i+2,3}=BPFamily(i).freqsLR(17,17);
    Temp{i+2,6}=BPFamily(i).nearfreqs(17,17);
    Temp{i+2,7}=BPFamily(i).nearfreqsLR(17,17);
end

xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Temp,'Summary','A1');

BPFamily=Reformat(BPFamily);

for i=1:12
   Label{1,1}=Names{i};
   if mod(i,2)==1
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'Freqs',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).freqs,'Freqs',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'RelFreqs',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).relfreqs,'RelFreqs',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'FreqsNear',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).nearfreqs,'FreqsNear',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'RelFreqsNear',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).nearrelfreqs,'RelFreqsNear',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'FreqsLR',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).freqsLR,'FreqsLR',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'FreqsLRNear',['A' num2str((i-1)/2*19+1) ':Q' num2str((i-1)/2*19+17)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).nearfreqsLR,'FreqsLRNear',['B' num2str((i-1)/2*19+2) ':R' num2str((i-1)/2*19+18)]);
   else
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'Freqs',['T' num2str((i-2)/2*19+1) ':AJ' num2str((i-2)/2*19+17)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],BPFamily(i).freqs,'Freqs',['U' num2str((i-2)/2*19+2) ':AK' num2str((i-2)/2*19+18)]);
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'RelFreqs',['T' num2stï¶2ÉŒé„BÃ@Hxe¦ƒÜ7º…ĞÃŒ;ŞeTtÒÇ ÖK`6<KU®åƒèTÎ~kæbÀgôm-i«ÎËË¢œ×Í®ûªiÃ#å9t›w½ Ó½İ°3#`J4¢ÖJMYÊ_V'Xü°…™«¬ê£şÖ‡PwZòûèåŞ&³ÅÃ(K=øÛ—Ê‹øäxEö®~Ù)ú©ÔN§ ÌÏgOš™µ…k48çğÀ«‡Ë-›\œ.¶	;ÀLÕãHã?g-»ñé\ÿd¥>Q˜9¬&ĞXÖúp’'Wı`	ÓÄ[äãÊXÉ¸ŠğÕıK(]¡ÏMo™¾	#ûõ±ŒX”vm,Ğ«JÓğş/s.¢Éñ¿ã¯±UØ
êQ	„íúõqÒö
Ë—“%ªâ`<­$®øA¦p@:ym—õdwdvjŒrŞƒm'+¼‰$õ25x½$9A“xOìvÜ<6~Ç {>_ç7?öú-Ìuû€#6(iÂcX€¸lÉ†0¯EÌ	óiè)­€³à2‘Íó5MµÖuäÂÙÄOÃt+òk>ßg‚‚DĞö„4,)ñ%/YdWÉ}D	 X´Û^K³`¹¨‹V¶¶ •¬}-¤³®¸MˆÌKtÆ¹Z|SË¼&@Ÿí’3™üï¶2ÉŒé„BÃ@Hxe¦ƒÜ7º…ĞÃŒ;ŞeTtÒÇ ÖK`6<KU®åƒèTÎ~kæbÀgôm-i«ÎËË¢œ×Í®ûªiÃ#å9t›w½ Ó½İ°3#`J4¢ÖJMYÊ_V'Xü°…™«¬ê£şÖ‡PwZòûèåŞ&³ÅÃ(K=øÛ—Ê‹øäxEö®~Ù)ú©ÔN§ ÌÏgOš™µ…k48çğÀ«‡Ë-›\œ.¶	;ÀLÕãHã?g-»ñé\ÿd¥>Q˜9¬&ĞXÖúp’'Wı`	ÓÄ[äãÊXÉ¸ŠğÕıK(]¡ÏMo™¾	#ûõ±ŒX”vm,Ğ«JÓğş/s.¢Éñ¿ã¯±UØ
êQ	„íúõqÒö
Ë—“%ªâ`<­$®øA¦p@:ym—õdwdvjŒrŞƒm'+¼‰$õ25x½$9A“xOìvÜ<6~Ç {>_ç7?öú-Ìuû€#6(iÂcX€¸lÉ†0¯EÌ	óiè)­€³à2‘Íó5MµÖuäÂÙÄOÃt+òk>ßg‚‚DĞö„4,)ñ%/YdWÉ}D	 X´Û^K³`¹¨‹V¶¶ •¬}-¤³®¸MˆÌKtÆ¹Z|SË¼&@Ÿí’3™ü…kç	»uu˜j(Ç­«Ûn±"
Z*FO ¨Nê9ÏlÃavrnø"ìcyF¶{ığßé}9j²ö*fâuÚF-!ÀS—º¤Í›ó¥›­rD²×ìhõ½¡¸öø<Åâî€‡ƒF N6Bf}32„h‰µHsé]”À¿IÊÈ‚3TäXè~„çuñ¡7b8w&ú#ú{i0”
eÅÉÏâ	Ò¯æ>eûM;ëm[†Âw–]è)tVo`Å b©±ú·ıÏ0fı7 ´jjåßMÈwDÛz#ïŞôé5R÷î$2ÍhoŒs&=í~íJ>ó¤ğôÈó\à&_Ì±“n„*¨fã¡#¬9’¢$½Y¨$o³[Å%İ²FQ))˜ˆ¼~4œR¾îÂh–÷²®t>©ˆë;š	c¿Ò	±åõ=wOWhàÛ8ó>ÁN‚áÃÏÜ¸.E› ”Ç}¬¼áâ®º»T7i´Œ19å&6a…ú‰y;Õé¯cM%ÂìÌ@óÑz­yÈIÀ/Ÿğq»³eAÑôRÂæ+ÀÍ>-Ğg9uğØ!ËßlySø6"bœÂhÔ{¯¯>ƒÂ;òt|%•Å«îf¾ĞFÿ{çÄ£ß˜¤)"s«Ñ‰1PHµ}{æ¯Øº‘Â…kç	»uu˜j(Ç­«Ûn±"
Z*FO ¨Nê9ÏlÃavrnø"ìcyF¶{ığßé}9j²ö*fâuÚF-!ÀS—º¤Í›ó¥›­rD²×ìhõ½¡¸öø<Åâî€‡ƒF N6Bf}32„h‰µHsé]”À¿IÊÈ‚3TäXè~„çuñ¡7b8w&ú#ú{i0”
eÅÉÏâ	Ò¯æ>eûM;ëm[†Âw–]è)tVo`Å b©±ú·ıÏ0fı7 ´jjåßMÈwDÛz#ïŞôé5R÷î$2ÍhoŒs&=í~íJ>ó¤ğôÈó\à&_Ì±“n„*¨fã¡#¬9’¢$½Y¨$o³[Å%İ²FQ))˜ˆ¼~4œR¾îÂh–÷²®t>©ˆë;š	c¿Ò	±åõ=wOWhàÛ8ó>ÁN‚áÃÏÜ¸.E› ”Ç}¬¼áâ®º»T7i´Œ19å&6a…ú‰y;Õé¯cM%ÂìÌ@óÑz­yÈIÀ/Ÿğq»³eAÑôRÂæ+ÀÍ>-Ğg9uğØ!ËßlySø6"bœÂhÔ{¯¯>ƒÂ;òt|%•Å«îf¾ĞFÿ{çÄ£ß˜¤)"s«Ñ‰1PHµ}{æ¯Øº‘ÂK#ÑZ º £ƒµš)ÎÊ ÙjÏŸ_HŠÏ:}Eb[=5ÃJ¯øäl‹
* vS¨ö5$‡o‡—ˆåÊè²÷A‡ôÂÜ’d°›+<jÛ[¹-øî‹xQÅ÷Kú†€û‹¯ÍëÇ¯B‹+Õ$·¾(V-ìØòjU¤š!E × S½SO¤A`äflŠahO!nYE‘Çµ?÷ûHMm©q’ë$[`û$–x&Aå—&U+–ÃÑõ9¹ƒªØÓ<‰sÜU›íŒkC¼õ,ù*İÁhˆŞ3hEÆîsªÆú½OB$`Åxß:Le˜(Ò5×„®“:ÃñXLËÊÉtùVÕü
Í G¡·÷(°¿Øw‰Ã§P Şà'4·s0‹HŸ7æ¤TA†İ5lsÁ½Ÿ<³—¶ƒ/gà£,Œ/d±ğP}½r‰ª
‰½µv÷¸M”ÒØb,¹'–ıUÜ6ÑşµB”i/Öğ—‚àe¸ª] |´(„wËOº4À÷"W”Íš˜Şşµ µsdÇ9Äà8vMw¢â)Ú”wWgQ³"CZÌ QtçTkÌ(5gÿáÖ~h9Î¸n ŸK†oKä.=ºWŠäÏöòÛE¬oĞ¾ï©pûé¥C2
—2›*sK#ÑZ º £ƒµš)ÎÊ ÙjÏŸ_HŠÏ:}Eb[=5ÃJ¯øäl‹
* vS¨ö5$‡o‡—ˆåÊè²÷A‡ôÂÜ’d°›+<jÛ[¹-øî‹xQÅ÷Kú†€û‹¯ÍëÇ¯B‹+Õ$·¾(V-ìØòjU¤š!E × S½SO¤A`äflŠahO!nYE‘Çµ?÷ûHMm©q’ë$[`û$–x&Aå—&U+–ÃÑõ9¹ƒªØÓ<‰sÜU›íŒkC¼õ,ù*İÁhˆŞ3hEÆîsªÆú½OB$`Åxß:Le˜(Ò5×„®“:ÃñXLËÊÉtùVÕü
Í G¡·÷(°¿Øw‰Ã§P Şà'4·s0‹HŸ7æ¤TA†İ5lsÁ½Ÿ<³—¶ƒ/gà£,Œ/d±ğP}½r‰ª
‰½µv÷¸M”ÒØb,¹'–ıUÜ6ÑşµB”i/Öğ—‚àe¸ª] |´(„wËOº4À÷"W”Íš˜Şşµ µsdÇ9Äà8vMw¢â)Ú”wWgQ³"CZÌ QtçTkÌ(5gÿáÖ~h9Î¸n ŸK†oKä.=ºWŠäÏöòÛE¬oĞ¾ï©pûé¥C2
—2›*sGZCëÜK²lÜ)“G½Ã™¤åQ²sY&ÛÔ-í³õ—ëMöó°>õ-w‚O\“÷Û,=?$Êİºø³Q"¥cŸå{ø†Y-W/:.Dğ[ñFeJ]E.›Æ­”}á×\Oêõ¼ÊÃ“y$—]$dt¤˜2 İÇ›!ÿ¼™´îQ¹ÿ¢1#(!©şR ÇCœXYn°Iq?«^ñ«­Ø›WÕ¡Äş÷ÚGFZ¥ÒKK².—Áó»²:´uJ˜dbÁ¬qù?xoA®uUÒíj#KkºÇ\Ş¢«)™ãÇwÿ,Î£J!¬Ê‹¤Aäü]0ª‘ FZ˜zK@›I9T6Ïé°!·×…Åfö»Û£©ÿo)ÅF[Î©êÊSEC3"1Q#Îúşı)	5ıÎ[˜¶8æ›»ÎÓjï8Íº5«ŒGÖŒŸk+z‡‚-}eÏ.O| æ…’O‘N’¨Ã=+êË8>&+pú‚õÀÂ*Xf?Éu­Â?ı:‚şL9=›‘„óÓtæ/M0g–LE5rZ›V{®ÄEæCiã¥V—p7]×ø¦ôX%fY‚nš„? ÒææØ\¨MÂÄ©8tÒw[™€dÆtB¡a $¼2ÓAîİBèaÆï2*:Gécë%0¥*×òAt*GZCëÜK²lÜ)“G½Ã™¤åQ²sY&ÛÔ-í³õ—ëMöó°>õ-w‚O\“÷Û,=?$Êİºø³Q"¥cŸå{ø†Y-W/:.Dğ[ñFeJ]E.›Æ­”}á×\Oêõ¼ÊÃ“y$—]$dt¤˜2 İÇ›!ÿ¼™´îQ¹ÿ¢1#(!©şR ÇCœXYn°Iq?«^ñ«­Ø›WÕ¡Äş÷ÚGFZ¥ÒKK².—Áó»²:´uJ˜dbÁ¬qù?xoA®uUÒíj#KkºÇ\Ş¢«)™ãÇwÿ,Î£J!¬Ê‹¤Aäü]0ª‘ FZ˜zK@›I9T6Ïé°!·×…Åfö»Û£©ÿo)ÅF[Î©êÊSEC3"1Q#Îúşı)	5ıÎ[˜¶8æ›»ÎÓjï8Íº5«ŒGÖŒŸk+z‡‚-}eÏ.O| æ…’O‘N’¨Ã=+êË8>&+pú‚õÀÂ*Xf?Éu­Â?ı:‚şL9=›‘„óÓtæ/M0g–LE5rZ›V{®ÄEæCiã¥V—p7]×ø¦ôX%fY‚nš„? ÒææØ\¨MÂÄ©8tÒw[™€dÆtB¡a $¼2ÓAîİBèaÆï2*:Gécë%0¥*×òAt*0%Qk¥…¦,å‚ƒ/«“,~ØBÌUVõQk„C¨;-ù}ôro“Y‰âaGÈ¥üíKåE|r¼¢{W¿ìıTj§SPæç³'ÍÌÚÂ5œsxàÕÃå–M.NÛƒŠ`¦êq¤ñŸ³–İøt®²RŸ(ÌVh,k}8É“«~°„iâÎ-òñe¬d\Eøj‡ş%”®Ğ‹ç¦‚·LßÈ‘}ˆz†XF,J»6èU¥ixÿ—9Ñd„ø_ñ×Ø*lõ¨Âvıú8i{…åËƒÉ’ÀÕq°?VWü S8 ¼¶Ëz²;2»‡5F9ïÁ¶ÇŞD’z™¼Ş’œ I¼'v;n¿c€=‰¿¯ó›{Hıæº}À”4á1,@\¶dC˜À×"æ„ù4ô”VÀYp™Èæùš¦Z
ë:rálâ§€aºù5Ÿï3AA"h{B–”ø’—,²«ä>¢,Úm¯¥Ù°\ÔE+[ÛJVƒ¾ÒÙ	WÜ¦Äæ¥?:ã\->‚©e^“ …ÏvÉ™LşÂÆ†µó„]ºº?L5”cÈV‚Õm·X-£' Ô§Fõœg¶á0G;97|ö±<#Û½€~øïô¾5Y{3ñ:m£–à©Ë]R‹æÍùÒÍÇV9"Ùkv´úŞ0%Qk¥…¦,å‚ƒ/«“,~ØBÌUVõQk„C¨;-ù}ôro“Y‰âaGÈ¥üíKåE|r¼¢{W¿ìıTj§SPæç³'ÍÌÚÂ5œsxàÕÃå–M.NÛƒŠ`¦êq¤ñŸ³–İøt®²RŸ(ÌVh,k}8É“«~°„iâÎ-òñe¬d\Eøj‡ş%”®Ğ‹ç¦‚·LßÈ‘}ˆz†XF,J»6èU¥ixÿ—9Ñd„ø_ñ×Ø*lõ¨Âvıú8i{…åËƒÉ’ÀÕq°?VWü S8 ¼¶Ëz²;2»‡5F9ïÁ¶ÇŞD’z™¼Ş’œ I¼'v;n¿c€=‰¿¯ó›{Hıæº}À”4á1,@\¶dC˜À×"æ„ù4ô”VÀYp™Èæùš¦Z
ë:rálâ§€aºù5Ÿï3AA"h{B–”ø’—,²«ä>¢,Úm¯¥Ù°\ÔE+[ÛJVƒ¾ÒÙ	WÜ¦Äæ¥?:ã\->‚©e^“ …ÏvÉ™LşÂÆ†µó„]ºº?L5”cÈV‚Õm·X-£' Ô§Fõœg¶á0G;97|ö±<#Û½€~øïô¾5Y{3ñ:m£–à©Ë]R‹æÍùÒÍÇV9"Ùkv´úŞP\{|bqwÀÃA# '›Æ	!3Ï¾™FB´ÄZ¤¹ô.Êàß$eäÁ*r,t?ÂóºÈøĞ1œ;}‰ı=Ç4J…²âägñéWsŸ²ı¦õ¶-Cá;ËÀ.ô‹:«7°bP±ÎÔ‰XıÛşg³şZ5µòï&ä;¢m½¿‘woúô©{w™f´7Æ9“v¿v%GŸŒyRxzäy.p“/f‹ØI7B•T³qĞÖIQ’Ş,T’·ÇÙ­â’ˆnHY£¨”LD^	?N)_wáO4Ë{Ù@W:ŸTÄõÍ„±_iÈ„Øòúƒ»§«@4pŒmœyŸ`§
Áp‡ágn\—¢MÊã>VŞpñ×İ]ª›4ZÆ˜œr›°BıÄ¼êô×±¦avf ùÇh½Ö<ä$à—Oø¸İÙ² hz)aóàÀfŸHè³œ:xìåo6‡¼©‰€@|1Na4ê½×WŸAa	y:¾’ÊâUw3_h£ÿ½sâÑoHLÒŠ¹ÕèÄ(¤Ú¾=óWlİHá´_ e“è™ÖÏx`kò¬Á‚p³‚ù‘*¡(Íw-vé49Œ‹¥*Ÿ™Ü(ø"„¥‘h- ]€ÑÁZÍgeĞlµçÏ/$ÅgÀ¾"±­ša¥Wür¶E…	P\{|bqwÀÃA# '›Æ	!3Ï¾™FB´ÄZ¤¹ô.Êàß$eäÁ*r,t?ÂóºÈøĞ1œ;}‰ı=Ç4J…²âägñéWsŸ²ı¦õ¶-Cá;ËÀ.ô‹:«7°bP±ÎÔ‰XıÛşg³şZ5µòï&ä;¢m½¿‘woúô©{w™f´7Æ9“v¿v%GŸŒyRxzäy.p“/f‹ØI7B•T³qĞÖIQ’Ş,T’·ÇÙ­â’ˆnHY£¨”LD^	?N)_wáO4Ë{Ù@W:ŸTÄõÍ„±_iÈ„Øòúƒ»§«@4pŒmœyŸ`§
Áp‡ágn\—¢MÊã>VŞpñ×İ]ª›4ZÆ˜œr›°BıÄ¼êô×±¦avf ùÇh½Ö<ä$à—Oø¸İÙ² hz)aóàÀfŸHè³œ:xìåo6‡¼©‰€@|1Na4ê½×WŸAa	y:¾’ÊâUw3_h£ÿ½sâÑoHLÒŠ¹ÕèÄ(¤Ú¾=óWlİHá´_ e“è™ÖÏx`kò¬Á‚p³‚ù‘*¡(Íw-vé49Œ‹¥*Ÿ™Ü(ø"„¥‘h- ]€ÑÁZÍgeĞlµçÏ/$ÅgÀ¾"±­ša¥Wür¶E…	¼¨âû%ıCÀıÅ×æuÏãW¡Å•j’[_«vlGyµ*RM‹"k©Ş©'Ò 0r36Å0´§·€¬¢ÈcŠÚŸûı@¤¦¶Ô8HÉu’-°}K<“ òK“ªËáèúœÜAUìiÆÄ9îªÍvFŒµ!Şz|•îŠ`4Dï´"c÷¹Ucı^'!	°b¼o¦2ÌéškB×Iáx,¦eådº|«j~…Œf‚ˆ£P„ÆÛ{XŠ_ì»ÄáS( oğšÛ9˜E¤ÏGsRª Ãî¶9àŞOÙKÛÁˆ—3ğQÆ²Xx¨¾^¹DU…ÄŞZ»ûÜ&Jil1–Ü“Ëşª n›hÿZ!Ê´køKAğ2\Õ.P>ZÂ»å']à{‘+ÊfMLo‡ZÚ€¹²ãbp»¦À;Qñ”?mÊ;ƒ«³¨Y‘!-æ?€(ºó?ª5f”š³ÿpk?´ç\‚7Ï%Ã‡·%r—ƒ‹İ+EòŒg{ùí"HÖ7hß÷T¸ıtÇÒ!…K™M•¹£­¡uî%Y¶î”I‚£‚ŞáLOÒò(Ù¹,“mê–öÙúËõ&ûyXŸú–;Á'®†Éûm–ån]üÙ(‘Ò±Ï¿ò=|Ã¬–«"ø­x£2¥®‡¢€†MãV¼¨âû%ıCÀıÅ×æuÏãW¡Å•j’[_«vlGyµ*RM‹"k©Ş©'Ò 0r36Å0´§·€¬¢ÈcŠÚŸûı@¤¦¶Ô8HÉu’-°}K<“ òK“ªËáèúœÜAUìiÆÄ9îªÍvFŒµ!Şz|•îŠ`4Dï´"c÷¹Ucı^'!	°b¼o¦2ÌéškB×Iáx,¦eådº|«j~…Œf‚ˆ£P„ÆÛ{XŠ_ì»ÄáS( oğšÛ9˜E¤ÏGsRª Ãî¶9àŞOÙKÛÁˆ—3ğQÆ²Xx¨¾^¹DU…ÄŞZ»ûÜ&Jil1–Ü“Ëşª n›hÿZ!Ê´køKAğ2\Õ.P>ZÂ»å']à{‘+ÊfMLo‡ZÚ€¹²ãbp»¦À;Qñ”?mÊ;ƒ«³¨Y‘!-æ?€(ºó?ª5f”š³ÿpk?´ç\‚7Ï%Ã‡·%r—ƒ‹İ+EòŒg{ùí"HÖ7hß÷T¸ıtÇÒ!…K™M•¹£­¡uî%Y¶î”I‚£‚ŞáLOÒò(Ù¹,“mê–öÙúËõ&ûyXŸú–;Á'®†Éûm–ån]üÙ(‘Ò±Ï¿ò=|Ã¬–«"ø­x£2¥®‡¢€†MãVÊ¾ğëŠ	®§õz^åáÉ<’Ë‰.²G:RLîãÍ^ŒLZwˆ¨ÜÑŒ”T)Ğã!N¬,7Ø¤¸ŸU¯øÕÖ@ìÍ+ˆêPbÿˆ{í##­RéÇ¥¥Y—Ëàƒù]YÚ:%ÈL2±`†GÖ¸ü¼· ×º*évµ‘¥5İc.oÑÕ”Ìñã»çÑ¥Ve‡EÒ rş.UÇH£N­GL½% Í¤*çtØÛëÂb3ûİíÑÔÿ·”b£-çTuå©¢¡Î‘…˜¨g}ÿş”„š~ç­L[óÍ]gƒé?µwœfİšUÆ#kÆÏµ½CÁ–¾²g—'>PóBÉ§È
'IÔáõe8}Áz`a,³ÀŸä¿ºVáŸ~ A¦œÍHÂùi:sŠ‡&˜3K¦¢9­M«=Wâ"ó!ƒ´ñR«K¸‹›®k|Sz¬³,A7ÍÂPissl.Ô&aHâTH:é»­LÀF2c:¡Ğ0^™é ÷n!ô0ãw£ô1ˆõ˜ÏRG•kù :•³D‡ßãš¹ğ}[KÚªóò²èçu³ë¾jÚ0ÂHyİæ]/ètcgo7ìÇÌ˜¨µÒBS–rÁÁ—ÕÉ?l¡@æ*«ú¨¿5Â!Ô–ü>z¹·É¬Dñ0Š#äÊ¾ğëŠ	®§õz^åáÉ<’Ë‰.²G:RLîãÍ^ŒLZwˆ¨ÜÑŒ”T)Ğã!N¬,7Ø¤¸ŸU¯øÕÖ@ìÍ+ˆêPbÿˆ{í##­RéÇ¥¥Y—Ëàƒù]YÚ:%ÈL2±`†GÖ¸ü¼· ×º*évµ‘¥5İc.oÑÕ”Ìñã»çÑ¥Ve‡EÒ rş.UÇH£N­GL½% Í¤*çtØÛëÂb3ûİíÑÔÿ·”b£-çTuå©¢¡Î‘…˜¨g}ÿş”„š~ç­L[óÍ]gƒé?µwœfİšUÆ#kÆÏµ½CÁ–¾²g—'>PóBÉ§È
'IÔáõe8}Áz`a,³ÀŸä¿ºVáŸ~ A¦œÍHÂùi:sŠ‡&˜3K¦¢9­M«=Wâ"ó!ƒ´ñR«K¸‹›®k|Sz¬³,A7ÍÂPissl.Ô&aHâTH:é»­LÀF2c:¡Ğ0^™é ÷n!ô0ãw£ô1ˆõ˜ÏRG•kù :•³D‡ßãš¹ğ}[KÚªóò²èçu³ë¾jÚ0ÂHyİæ]/ètcgo7ìÇÌ˜¨µÒBS–rÁÁ—ÕÉ?l¡@æ*«ú¨¿5Â!Ô–ü>z¹·É¬Dñ0Š#ä³ßİMıK)6ÚrNUW*êœYˆ‰qÆĞ÷ïOI¨éwŞêÀ´Å1ßÜu6˜şS{ÇiÖ­Ye<²fü\[Ñë8lé+{vyâ5/”|Š¬p’DîYQ_Æñ1áXÓ¬VÁ2üIş«kşéçĞyôgÊéÙŒ$œŸ¦3§xqh‚9³d*ª‘ÓÚ´Ús%.22H/µº„»¸éºÆ7¥Ç*1ËtÓ ü•67ÇæBm†$NÅ¤“¾ÛÊl$3¦
!á•™rßèB3îx—QÑ9JƒX/Ùğ,uT¹–¢S9Ktø=®™‹ŸÑ·µ¤­://‹>p^7»î«¦#Œ”çĞmŞõ‚N7vövÃ~ÌŒ€)ÑˆZ+-4e)|Y,`ñÃ
d®²ªú[#BİiÉï£—{›ÌJ£8B.õào_*/â“ã-Ø»úe§è§R;‚2?Ÿ=ifÖ®ÑàœÃ¯.·lrqºØ$Pì 3U#ÿœµìÆ§sı“•úDaæ°š@cYëÃ	H\õƒ%Lwn‘(c%ã*ÂW;ô/¡t…^<7¼eú&@ìCÔ3Ä2bQÚµ±@¯*MÃû¿Ì¹ˆ&#ÄÿzŒ¿Æ^Ta+¨G%¶ë×ÇI{Ø+,_L–.¨.ˆƒıñ´’¸rà™Âyèä³ßİMıK)6ÚrNUW*êœYˆ‰qÆĞ÷ïOI¨éwŞêÀ´Å1ßÜu6˜şS{ÇiÖ­Ye<²fü\[Ñë8lé+{vyâ5/”|Š¬p’DîYQ_Æñ1áXÓ¬VÁ2üIş«kşéçĞyôgÊéÙŒ$œŸ¦3§xqh‚9³d*ª‘ÓÚ´Ús%.22H/µº„»¸éºÆ7¥Ç*1ËtÓ ü•67ÇæBm†$NÅ¤“¾ÛÊl$3¦
!á•™rßèB3îx—QÑ9JƒX/Ùğ,uT¹–¢S9Ktø=®™‹ŸÑ·µ¤­://‹>p^7»î«¦#Œ”çĞmŞõ‚N7vövÃ~ÌŒ€)ÑˆZ+-4e)|Y,`ñÃ
d®²ªú[#BİiÉï£—{›ÌJ£8B.õào_*/â“ã-Ø»úe§è§R;‚2?Ÿ=ifÖ®ÑàœÃ¯.·lrqºØ$Pì 3U#ÿœµìÆ§sı“•úDaæ°š@cYëÃ	H\õƒ%Lwn‘(c%ã*ÂW;ô/¡t…^<7¼eú&@ìCÔ3Ä2bQÚµ±@¯*MÃû¿Ì¹ˆ&#ÄÿzŒ¿Æ^Ta+¨G%¶ë×ÇI{Ø+,_L–.¨.ˆƒıñ´’¸rà™Âyèäµ]Ö“İ‘Ù=4¨1Êy¶8®ğ&’ÔËÔàõnäMâ=±ÛqóØøìùHü}ßüØCê·0×íØ  ¤	aâ²%:À¾~1'Ì§¡§´Î‚ËD6Ï×4ÕRX×‘g?Ó­0È¯ù|Ÿ	
AÛÒ°¤Ä—¼d‘]%ô%€`Ñncx-Í€å¢.ZÙÚf€T²ôµÎN¸â6= f0/ıÑçjñL-óšt -|´KÎdò66¬'ì
ÔÕıab¨¡C¶¬n»ÅŠ(h©= >85ªç<³‡9ÚÉ¹á‹°åÙîôÃ§÷å¨ÉÚ«˜‰×iµ„ O]fè’Z4oÎ—n>¶ÊÉ^³£Õ÷†âÚãó‹» 8Ùt4N™yöÍ4Ê¢%Ö"Í¥wQ~ ÿ&)#ÎP‘c¡ûwÔ@Æ‡ŞˆáÜ™èKŒèï9¦ÁP*”'?;ˆ'H¿šû”í7í¬·m
ßYv¡_¤ĞY½ƒŠu¦NÄêßö?Ã˜õß€Ğª!¨•7!ßmëı¼{Ó§×Hİ»“È4£½1nÌA˜ô´ûµ+9údÌ“ÂÓ#Ïs›|1[ÄNj¸ª š{„°æHŠ’ôf¡’¼<În—DtCÊE¥¤`"òJøÑpJùº/x¢YŞËºÒù¤"®ïµ]Ö“İ‘Ù=4¨1Êy¶8®ğ&’ÔËÔàõnäMâ=±ÛqóØøìùHü}ßüØCê·0×íØ  ¤	aâ²%:À¾~1'Ì§¡§´Î‚ËD6Ï×4ÕRX×‘g?Ó­0È¯ù|Ÿ	
AÛÒ°¤Ä—¼d‘]%ô%€`Ñncx-Í€å¢.ZÙÚf€T²ôµÎN¸â6= f0/ıÑçjñL-óšt -|´KÎdò66¬'ì
ÔÕıab¨¡C¶¬n»ÅŠ(h©= >85ªç<³‡9ÚÉ¹á‹°åÙîôÃ§÷å¨ÉÚ«˜‰×iµ„ O]fè’Z4oÎ—n>¶ÊÉ^³£Õ÷†âÚãó‹» 8Ùt4N™yöÍ4Ê¢%Ö"Í¥wQ~ ÿ&)#ÎP‘c¡ûwÔ@Æ‡ŞˆáÜ™èKŒèï9¦ÁP*”'?;ˆ'H¿šû”í7í¬·m
ßYv¡_¤ĞY½ƒŠu¦NÄêßö?Ã˜õß€Ğª!¨•7!ßmëı¼{Ó§×Hİ»“È4£½1nÌA˜ô´ûµ+9údÌ“ÂÓ#Ïs›|1[ÄNj¸ª š{„°æHŠ’ôf¡’¼<În—DtCÊE¥¤`"òJøÑpJùº/x¢YŞËºÒù¤"®ïèîRİ¤Ñ2Æä”›Ø„ê'æíT§¿5•³3Í?Fëµæ!'¿|ÂÇíÎ–EÓK	›¯ 6û´@BŸåÔÁ3`‡,³9äMMâÛˆˆq
£Qï½¾ú
{ìHXÈÓñ•T¯º›ùBıï~C@b2¦ˆPÌ­F'Æ@!Õöí™¿bëFz§ı(»˜DÏ´~Æk [“g„›œÈT	í@i¾k±K§Éa\,UùÌäFÁÿ!,Dk0èŒÖj¦8+ƒf«u<~!9(>ëö‰mõÔ+½â_³-*L¨€ÚM¡Ú×`v¾^"”+£ËŞuÒSpK’Án®`ğ¨moå¶à[<^¸/.àEß/éî/¾6¬{¿
-