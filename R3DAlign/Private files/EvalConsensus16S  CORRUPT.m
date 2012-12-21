% -------------------------------------- Load Jesse/Ryan Composite 3D to 3D
% alignment usign Ryan's better alignment (seedAcc=25);
a = 1;
Al(a).Name           = 'Composite2';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load Jesse/Ryan Composite 3D to 3D
% alignment using faster alignment
a = a+1;
Al(a).Name           = 'Composite1';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load Ryan's 3D to 3D alignment
a = a+1;
Al(a).Name           = 'Ryan';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load Ryan's 3D to 3D alignment
a = a+1;
Al(a).Name           = 'Ryan2';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load Jesse's 3D to 3D alignment
a = a+1;
Al(a).Name           = 'Jesse';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load SARA 3D to 3D alignment
a = a+1;
Al(a).Name           = 'SARA';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load Ryan/SARA 3D to 3D alignment
a = a+1;
Al(a).Name           = 'Ryan/SARA';
Al(a).SheetName='RyanSara';

% -------------------------------------- Load ARTS 3D to 3D alignment
a = a+1;
Al(a).Name           = 'ARTS';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load DIAL 3D to 3D alignment
a = a+1;
Al(a).Name           = 'DIAL';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load SARSA 3D to 3D alignment
a = a+1;
Al(a).Name           = 'SARSA';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load NW alignment
a = a+1;
Al(a).Name           = 'NW';
Al(a).SheetName=Al(a).Name;

% -------------------------------------- Load FoldAlign alignment
a = a+1;
Al(a).Name           = 'FoldAlign';
Al(a).SheetName='FoldAl';

T=cell(17,length(Al)+1);

for k = 1:length(Al)
    clear a;
    clear b;
    Al(k).SheetName
    [a b]=xlsread([pwd filesep 'rAlignmentBasePairComparison9 1J5E 2AVY'],Al(k).SheetName);
   for i=1:length(b(:,2))
      b(i,2)=lower(strtrim(b(i,2)));
      b(i,5)=lower(strtrim(b(i,5)));
   end
   ct=0;
   for i=1:length(b(:,2))
      if isequal(b{i,2},'tsh')
         b{i,2}='ths';
      elseif isequal(b{i,2},'csh')
         b{i,2}='chs';
      elseif isequal(b{i,2},'tsw')
         b{i,2}='tws';
      elseif isequal(b{i,2},'csw')
         b{i,2}='cws';
      elseif isequal(b{i,2},'thw')
         b{i,2}='twh';
      elseif isequal(b{i,2},'chw')
         b{i,2}='cwh';
      elseif isequal(b{i,2},'ntsh')
         b{i,2}='nths';
      elseif isequal(b{i,2},'ncsh')
         b{i,2}='nchs';
      elseif isequal(b{i,2},'ntsw')
         b{i,2}='ntws';
      elseif isequal(b{i,2},'ncsw')
         b{i,2}='ncws';
      elseif isequal(b{i,2},'nthw')
         b{i,2}='ntwh';
      elseif isequal(b{i,2},'nchw')
         b{i,2}='ncwh';
      elseif ~isempty(b{i,2}) && isequal(b{i,2}(1),'s')
         b{i,2}='';
      end
   end
   for i=1:length(b(:,5))
      if isequal(b{i,5},'tsh')
         b{i,5}='ths';
      elseif isequal(b{i,5},'csh')
         b{i,5}='chs';
      elseif isequal(b{i,5},'tsw')
         b{i,5}='tws';
      elseif isequal(b{i,5},'csw')
         b{i,5}='cws';
      elseif isequal(b{i,5},'thw')
         b{i,5}='twh';
      elseif isequal(b{i,5},'chw')
         b{i,5}='cwh';
      elseif isequal(b{i,5},'ntsh')
         b{i,5}='nths';
      elseif isequal(b{i,5},'ncsh')
         b{i,5}='nchs';
      elseif isequal(b{i,5},'ntsw')
         b{i,5}='ntws';
      elseif isequal(b{i,5},'ncsw')
         b{i,5}='ncws';
      elseif isequal(b{i,5},'nthw')
         b{i,5}='ntwh';
      elseif isequal(b{i,5},'nchw')
         b{i,5}='ncwh';
      elseif ~isempty(b{i,5}) && isequal(b{i,5}(1),'s')
         b{i,5}='';
      end
   end

   numbp1=0;
   numnbp1=0;
   numnothing1=0;
   for i=1:length(b(:,2))
      len=length(b{i,2});
      if len==3
         numbp1=numbp1+1;
      elseif len==4 && ~isequal(b{i,2},'----')
         numnbp1=numnbp1+1;
      else 
         numnothing1=numnothing1+1;
      end
   end
%    [numbp1 numnbp1 numnothing1]

   numbp2=0;
   numnbp2=0;
   numnothing2=0;
   for i=1:length(b(:,5))
      len=length(b{i,5});
      if len==3
         numbp2=numbp2+1;
      elseif len==4 && ~isequal(b{i,5},'----')
         numnbp2=numnbp2+1;
      else 
         numnothing2=numnothing2+1;
      end
   end
%   [numbp2 numnbp2 numnothing2]

%We're comparing columns two and five.  
% Different Scenarios:
% 1) same basepairs
% 2) Different basepairs
% 3) Basepair in A, correct near bp in B
% 4) Bp in A, different near bp in B
% 5) Bp in A, aligned with 2 nt's making no bp in B
% 6) Bp in A, neither is aligned
% 7) Bp in A, only one is aligned
% 8) Basepair in B, correct near bp in A
% 9) Bp in B, different near bp in B
% 10) Bp in B, aligned with 2 nt's making no bp in A 
% 11) Bp in B, neither is aligned
% 12) Bp in B, only one is aligned
% 13) Single nucleotide in A aligned with nothing
% 14) Single nucleotide in B aligned with nothing
% 15) Single nucleotides aligned
% 16) Nbp in A with Nbp in B
% 17) Nbp in A, not with Nbp in B
% 18) Nbp in B, not with Nbp in A
   samebp = 0;
   diffbp = 0;
   inAnotB = 0;
   AwithNear = 0;
   AwithWrongNear = 0;
   AwithNoNT = 0;
   AwithOneNT = 0;
   inBnotA = 0;
   BwithNear = 0;
   BwithWrongNear = 0;
   BwithNoNT = 0;
   BwithOneNT = 0;
   SingleAwithNot = 0;
   SingleBwithNot = 0;
   SinglesAligned = 0;
   twoNears = 0;
   NearInAnotWithNearinB = 0;
   NearInBnotWithNearinA = 0;
   for i=1:length(b(:,2))
      if length(b{i,2})==3    
         if length(b{i,5})==3
            if isequal(b{i,2},b{i,5})
               samebp=samebp+1;                           % Case 1
            else
               diffbp=diffbp+1;                           % Case 2
            end
         elseif length(b{i,5})==4
            if isequal(b{i,5},'----')
               inAnotB=inAnotB+1;                        % Case 5
            elseif isequal(strcat('n',b{i,2}),b{i,5}) 
               AwithNear=AwithNear+1;                    % Case 3
            else
               AwithWrongNear=AwithWrongNear+1;          % Case 4
            end
         elseif isequal(b{i,5},'')
            if isequal(b{i,4},'---') && isequal(b{i,6},'---')
               AwithNoNT = AwithNoNT + 1;                % Case 6
            elseif isequal(b{i,4},'---') || isequal(b{i,6},'---')
               AwithOneNT = AwithOneNT + 1;              % Case 7
            end
         end
      elseif length(b{i,5})==3
         if length(b{i,2})==4
            if isequal(b{i,2},'----')
               inBnotA=inBnotA+1;                        % Case 10
            elseif isequal(strcat('n',b{i,5}),b{i,2})
               BwithNear=BwithNear+1;                    % Case 8
            else
               BwithWrongNear=BwithWrongNear+1;          % Case 9
            end
         elseif isequal(b{i,2},'')
            if isequal(b{i,1},'---') && isequal(b{i,3},'---')
               BwithNoNT = BwithNoNT + 1;                % Case 11
            elseif isequal(b{i,1},'---') || isequal(b{i,3},'---')
               BwithOneNT = BwithOneNT + 1;              % Case 12
            end
         end
      elseif isequal(b{i,4},'---') && isequal(b{i,6},'---')
         SingleAwithNot = SingleAwithNot + 1;            % Case 13
      elseif isequal(b{i,1},'---') && isequal(b{i,3},'---')
         SingleBwithNot = SingleBwithNot + 1;            % Case 14
      elseif isequal(b{i,3},'---') && isequal(b{i,6},'---') && ~isequal(b{i,1},'---') && ~isequal(b{i,4},'---')
         SinglesAligned = SinglesAligned + 1;            % Case 15
      elseif ~isempty(b{i,2}) && isequal(b{i,2}(1),'n')
         if isempty(b{i,5})     
             NearInAnotWithNearinB = NearInAnotWithNearinB + 1;
         elseif isequal(b{i,5}(1),'n')
             twoNears = twoNears+1;
         else
             NearInAnotWithNearinB = NearInAnotWithNearinB + 1;
         end
      i
      elseif isequal(b{i,5}(Ì²İ4@¥ÍÍ±¹P›„!‰Sq é¤ï¶2ÉŒé„BÃ@Hxe¦ƒÜ7º…ĞÃŒ;ŞeTtÒÇ ÖK`6<KU®åƒèTÎ~kæbÀgôm-i«ÎËË¢œ×Í®ûªiÃ#å9t›w½ Ó½İ°3#`J4¢ÖJMYÊ_V'Xü°…™«¬ê£şÖ‡PwZòûèåŞ&³ÅÃ(K=øÛ—Ê‹øäxEö®~Ù)ú©ÔN§ ÌÏgOš™µ…k48çğÀ«‡Ë-›\œ.¶	;ÀLÕãHã?g-»ñé\ÿd¥>Q˜9¬&ĞXÖúp’'Wı`	ÓÄ[äãÊXÉ¸ŠğÕıK(]¡ÏMo™¾	#ûõ±ŒX”vm,Ğ«JÓğş/s.¢Éñ¿ã¯±UØ
êQ	„íúõqÒö
Ë—“%ªâ`<­$®øA¦p@:ym—õdwdvjŒrŞƒm'+¼‰$õ25x½$9A“xOìvÜ<6~Ç {>_ç7?öú-Ìuû€#6(iÂcX€¸lÉ†0¯EÌ	óiè)­€³à2‘Íó5MµÖuäÂÙÄOÃt+òk>ßg‚‚DĞö„4,)ñ%/YdWÉ}D	 X´Û^K³`¹¨‹V¶¶ •¬}-¤³®¸MˆÌ²İ4@¥ÍÍ±¹P›„!‰Sq é¤ï¶2ÉŒé„BÃ@Hxe¦ƒÜ7º…ĞÃŒ;ŞeTtÒÇ ÖK`6<KU®åƒèTÎ~kæbÀgôm-i«ÎËË¢œ×Í®ûªiÃ#å9t›w½ Ó½İ°3#`J4¢ÖJMYÊ_V'Xü°…™«¬ê£şÖ‡PwZòûèåŞ&³ÅÃ(K=øÛ—Ê‹øäxEö®~Ù)ú©ÔN§ ÌÏgOš™µ…k48çğÀ«‡Ë-›\œ.¶	;ÀLÕãHã?g-»ñé\ÿd¥>Q˜9¬&ĞXÖúp’'Wı`	ÓÄ[äãÊXÉ¸ŠğÕıK(]¡ÏMo™¾	#ûõ±ŒX”vm,Ğ«JÓğş/s.¢Éñ¿ã¯±UØ
êQ	„íúõqÒö
Ë—“%ªâ`<­$®øA¦p@:ym—õdwdvjŒrŞƒm'+¼‰$õ25x½$9A“xOìvÜ<6~Ç {>_ç7?öú-Ìuû€#6(iÂcX€¸lÉ†0¯EÌ	óiè)­€³à2‘Íó5MµÖuäÂÙÄOÃt+òk>ßg‚‚DĞö„4,)ñ%/YdWÉ}D	 X´Û^K³`¹¨‹V¶¶ •¬}-¤³®¸MˆÌKtÆ¹Z|SË¼&@Ÿí’3™ü…kç	»uu˜j(Ç­«Ûn±"
Z*FO ¨Nê9ÏlÃavrnø"ìcyF¶{ığßé}9j²ö*fâuÚF-!ÀS—º¤Í›ó¥›­rD²×ìhõ½¡¸öø<Åâî€‡ƒF N6Bf}32„h‰µHsé]”À¿IÊÈ‚3TäXè~„çuñ¡7b8w&ú#ú{i0”
eÅÉÏâ	Ò¯æ>eûM;ëm[†Âw–]è)tVo`Å b©±ú·ıÏ0fı7 ´jjåßMÈwDÛz#ïŞôé5R÷î$2ÍhoŒs&=í~íJ>ó¤ğôÈó\à&_Ì±“n„*¨fã¡#¬9’¢$½Y¨$o³[Å%İ²FQ))˜ˆ¼~4œR¾îÂh–÷²®t>©ˆë;š	c¿Ò	±åõ=wOWhàÛ8ó>ÁN‚áÃÏÜ¸.E› ”Ç}¬¼áâ®º»T7i´Œ19å&6a…ú‰y;Õé¯cM%ÂìÌ@óÑz­yÈIÀ/Ÿğq»³eAÑôRÂæ+ÀÍ>-Ğg9uğØ!ËßlySø6"bœÂhÔ{¯¯>ƒÂ;òt|%•Å«îf¾ĞFÿ{çÄ£ßÌKtÆ¹Z|SË¼&@Ÿí’3™ü…kç	»uu˜j(Ç­«Ûn±"
Z*FO ¨Nê9ÏlÃavrnø"ìcyF¶{ığßé}9j²ö*fâuÚF-!ÀS—º¤Í›ó¥›­rD²×ìhõ½¡¸öø<Åâî€‡ƒF N6Bf}32„h‰µHsé]”À¿IÊÈ‚3TäXè~„çuñ¡7b8w&ú#ú{i0”
eÅÉÏâ	Ò¯æ>eûM;ëm[†Âw–]è)tVo`Å b©±ú·ıÏ0fı7 ´jjåßMÈwDÛz#ïŞôé5R÷î$2ÍhoŒs&=í~íJ>ó¤ğôÈó\à&_Ì±“n„*¨fã¡#¬9’¢$½Y¨$o³[Å%İ²FQ))˜ˆ¼~4œR¾îÂh–÷²®t>©ˆë;š	c¿Ò	±åõ=wOWhàÛ8ó>ÁN‚áÃÏÜ¸.E› ”Ç}¬¼áâ®º»T7i´Œ19å&6a…ú‰y;Õé¯cM%ÂìÌ@óÑz­yÈIÀ/Ÿğq»³eAÑôRÂæ+ÀÍ>-Ğg9uğØ!ËßlySø6"bœÂhÔ{¯¯>ƒÂ;òt|%•Å«îf¾ĞFÿ{çÄ£ßò#UB;PšïZìÒirKU>3¹QğÿDK#ÑZ º £ƒµš)ÎÊ ÙjÏŸ_HŠÏ:}Eb[=5ÃJ¯øäl‹
* vS¨ö5$‡o‡—ˆåÊè²÷A‡ôÂÜ’d°›+<jÛ[¹-øî‹xQÅ÷Kú†€û‹¯ÍëÇ¯B‹+Õ$·¾(V-ìØòjU¤š!E × S½SO¤A`äflŠahO!nYE‘Çµ?÷ûHMm©q’ë$[`û$–x&Aå—&U+–ÃÑõ9¹ƒªØÓ<‰sÜU›íŒkC¼õ,ù*İÁhˆŞ3hEÆîsªÆú½OB$`Åxß:Le˜(Ò5×„®“:ÃñXLËÊÉtùVÕü
Í G¡·÷(°¿Øw‰Ã§P Şà'4·s0‹HŸ7æ¤TA†İ5lsÁ½Ÿ<³—¶ƒ/gà£,Œ/d±ğP}½r‰ª
‰½µv÷¸M”ÒØb,¹'–ıUÜ6ÑşµB”i/Öğ—‚àe¸ª] |´(„wËOº4À÷"W”Íš˜Şşµ µsdÇ9Äà8vMw¢â)Ú”wWgQ³"CZÌ QtçTkÌ(5gÿáÖ~h9Î¸n ŸK†oKä.=ºWŠäÏò#UB;PšïZìÒirKU>3¹QğÿDK#ÑZ º £ƒµš)ÎÊ ÙjÏŸ_HŠÏ:}Eb[=5ÃJ¯øäl‹
* vS¨ö5$‡o‡—ˆåÊè²÷A‡ôÂÜ’d°›+<jÛ[¹-ø