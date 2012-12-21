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
      xlswrite([File1.Filename '_' File2.Filename '_alignment_data.xlsx'],Label,'RelFreqs',['T' num2st�2Ɍ�B�@Hxe���7���Ì;�eTt��� �K`6<KU���T�~�k�b�g�m-i���ˢ��ͮ��i�#�9t�w��Ӎ��ݰ3#`J4��JMY�_V'X����������PwZ�����&���(��K=�ۗʋ��xE��~�)���N����gO����k48������-�\�.�	;�L��H�?g-���\�d�>Q�9�&�X��p�'W�`	�ĝ[���Xɸ����K(]��Mo��	�#����X�vm,ЫJ���/s.���㯱U�
�Q	����q��
˗�%���`<�$��A�p@:ym��dwdvj�rރm'�+��$�25x�$9A�xO�v�<6~� {>_�7?���-�u��#6(i�cX��lɆ0��E�	�i�)����2���5M��u����O�t+�k>�g��D���4,)�%/YdW�}D	 X��^K�`���V�� ��}-����M��KtƹZ|S˼&@��3���2Ɍ�B�@Hxe���7���Ì;�eTt��� �K`6<KU���T�~�k�b�g�m-i���ˢ��ͮ��i�#�9t�w��Ӎ��ݰ3#`J4��JMY�_V'X����������PwZ�����&���(��K=�ۗʋ��xE��~�)���N����gO����k48������-�\�.�	;�L��H�?g-���\�d�>Q�9�&�X��p�'W�`	�ĝ[���Xɸ����K(]��Mo��	�#����X�vm,ЫJ���/s.���㯱U�
�Q	����q��
˗�%���`<�$��A�p@:ym��dwdvj�rރm'�+��$�25x�$9A�xO�v�<6~� {>_�7?���-�u��#6(i�cX��lɆ0��E�	�i�)����2���5M��u����O�t+�k>�g��D���4,)�%/YdW�}D	 X��^K�`���V�� ��}-����M��KtƹZ|S˼&@��3����k�	�uu�j(ǐ���n�"
Z*FO �N��9�l�a�vrn�"�cyF�{����}9j��*f�u�F-!�S���͛󥛏�rD���h������<���F N6�Bf�}3�2�h��Hs�]���I���3T�X�~��u��7b8w&�#�{�i0�
e����	ү�>e�M;�m[��w��]�)tVo`Šb�������0f�7 �jj��M�wD�z#����5R��$2�ho�s&=�~�J�>�����\�&_���n�*�f��#�9��$�Y�$o��[�%ݐ�FQ))���~4�R����h�����t>���;�	c�Ґ	���=wOW�h��8�>�N����ܸ.E� ��}�������T7i��19�&6a���y;��cM%���@��z�y�I�/��q��eA��R��+���>-��g9u��!��lyS��6"b��h�{��>��;�t|%�ū�f��F�{�ģ����)"s�щ1PH�}{�غ��k�	�uu�j(ǐ���n�"
Z*FO �N��9�l�a�vrn�"�cyF�{����}9j��*f�u�F-!�S���͛󥛏�rD���h������<���F N6�Bf�}3�2�h��Hs�]���I���3T�X�~��u��7b8w&�#�{�i0�
e����	ү�>e�M;�m[��w��]�)tVo`Šb�������0f�7 �jj��M�wD�z#����5R��$2�ho�s&=�~�J�>�����\�&_���n�*�f��#�9��$�Y�$o��[�%ݐ�FQ))���~4�R����h�����t>���;�	c�Ґ	���=wOW�h��8�>�N����ܸ.E� ��}�������T7i��19�&6a���y;��cM%���@��z�y�I�/��q��eA��R��+���>-��g9u��!��lyS��6"b��h�{��>��;�t|%�ū�f��F�{�ģ����)"s�щ1PH�}{�غ��K#�Z � ����)�ʠ�jϟ_H��:�}Eb[=5�J���l�
*�vS��5$��o�������A���ܒd��+<j�[�-���xQ��K��������ǯB�+�$��(V-�؎�jU��!E � S�SO�A`�fl�ahO!nYE���?���HMm�q���$[`�$�x&A�&U+����9�����<��s�U��kC��,�*��h��3hE��s����OB$`�x�:Le�(�5ׄ��:��XL���t�V��
� G����(���w�çP ��'4�s0�H�7��TA��5ls���<����/g�,�/d��P}�r��
���v��M���b,�'��U�6���B�i/���e��]�|�(�w�O�4��"W�͚���� �sd�9��8vM�w��)ڔwWgQ�"CZ� Qt�Tk�(5g���~h9��n �K�oK�.=�W������E��oо�p�鎥C2
�2�*sK#�Z � ����)�ʠ�jϟ_H��:�}Eb[=5�J���l�
*�vS��5$��o�������A���ܒd��+<j�[�-���xQ��K��������ǯB�+�$��(V-�؎�jU��!E � S�SO�A`�fl�ahO!nYE���?���HMm�q���$[`�$�x&A�&U+����9�����<��s�U��kC��,�*��h��3hE��s����OB$`�x�:Le�(�5ׄ��:��XL���t�V��
� G����(���w�çP ��'4�s0�H�7��TA��5ls���<����/g�,�/d��P}�r��
���v��M���b,�'��U�6���B�i/���e��]�|�(�w�O�4��"W�͚���� �sd�9��8vM�w��)ڔwWgQ�"CZ� Qt�Tk�(5g���~h9��n �K�oK�.=�W������E��oо�p�鎥C2
�2�*sGZC��K�l�)�G�Ù���Q�sY&��-����M��>�-w�O\���,=?$�ݺ��Q"�c��{��Y-W/:.D�[�FeJ]E.�ƭ�}��\O����Óy$�]$d�t��2 �Ǜ!�����Q���1#(!��R��C�XYn�Iq?�^����؛Wա����GFZ�ҏKK�.���:�uJ��db���q�?xoA�uU��j#Kk��\ޢ�)���w�,ΣJ!����A��]0��� F�Z��zK@�I9T6��!�ׅ�f��ۣ��o)�F[Ω��SEC�3"1Q#����)	5��[��8曻��j�8ͺ5��G֌�k+z��-}e�.O|�慒O�N���=+��8>&+p�����*Xf�?�u��?�:��L9=�����t�/M0g�LE5rZ�V{��E�Ci�V�p7]����X%fY�n��?�����\�Mĩ8�t�w[���d�tB�a $�2�A��B�a��2*:G�c�%0���*��At*GZC��K�l�)�G�Ù���Q�sY&��-����M��>�-w�O\���,=?$�ݺ��Q"�c��{��Y-W/:.D�[�FeJ]E.�ƭ�}��\O����Óy$�]$d�t��2 �Ǜ!�����Q���1#(!��R��C�XYn�Iq?�^����؛Wա����GFZ�ҏKK�.���:�uJ��db���q�?xoA�uU��j#Kk��\ޢ�)���w�,ΣJ!����A��]0��� F�Z��zK@�I9T6��!�ׅ�f��ۣ��o)�F[Ω��SEC�3"1Q#����)	5��[��8曻��j�8ͺ5��G֌�k+z��-}e�.O|�慒O�N���=+��8>&+p�����*Xf�?�u��?�:��L9=�����t�/M0g�LE5rZ�V{��E�Ci�V�p7]����X%fY�n��?�����\�Mĩ8�t�w[���d�tB�a $�2�A��B�a��2*:G�c�%0���*��At*0%Qk���,傃/��,~�B��UV�Qk�C�;-�}�ro�Y��aGȥ��K�E|r��{W���Tj�SP��'����5�sx����M.Nۃ�`��q�񟳖��t��R�(�Vh,k}8ɓ�~��i��-��e�d\E�j��%��Ћ禂�L�ȑ}�z�XF,J�6�U�ix��9�d��_�����*l���v��8i{��˃ɒ��q�?�VW� S8 ����z�;2��5F9�����D�z������I�'v;n�c�=����{H��}��4�1,@\�dC��׏"��4��V�Yp������Z
�:r�l⧀a��5��3AA"h{B�����,���>�,�m����\�E+[��JV����	Wܦ��?:�\->��e^���ρvəL��Ɔ��]���?L5�c�V��m�X-�' ��F��g��0G;97|��<#۽�~����5Y{3�:m����]R�������V9"�kv���0%Qk���,傃/��,~�B��UV�Qk�C�;-�}�ro�Y��aGȥ��K�E|r��{W���Tj�SP��'����5�sx����M.Nۃ�`��q�񟳖��t��R�(�Vh,k}8ɓ�~��i��-��e�d\E�j��%��Ћ禂�L�ȑ}�z�XF,J�6�U�ix��9�d��_�����*l���v��8i{��˃ɒ��q�?�VW� S8 ����z�;2��5F9�����D�z������I�'v;n�c�=����{H��}��4�1,@\�dC��׏"��4��V�Yp������Z
�:r�l⧀a��5��3AA"h{B�����,���>�,�m����\�E+[��JV����	Wܦ��?:�\->��e^���ρvəL��Ɔ��]���?L5�c�V��m�X-�' ��F��g��0G;97|��<#۽�~����5Y{3�:m����]R�������V9"�kv���P\{|�bqw��A# '���	!3Ͼ�FB��Z���.���$e��*r,t?�����1�;}��=�4J����g��Ws�������-C�;��.�:�7�bP��ԉX���g��Z5���&�;�m���wo���{w�f�7ƍ9��v�v%G��yRxz�y.p�/f��I7B�T�q���IQ��,T����٭⒈nHY���LD^	?N)_w�O4�{�@W:�T���̈́�_iȄ��������@4p�m�y�`��
�p��gn\��M��>V�p���]��4ZƘ�r��B�ļ���ױ�avf���h��<�$��O���ٲ�hz)a���f�H賜:x��o6�����@|1Na4��W�Aa�	y:����Uw3_h���s��oHL������(�ھ=�Wl�H�_ e����x`k���p����*�(�w-v�49���*���(�"���h- ]���Z�ge�l����/$�g���"����a�W�r�E�	P\{|�bqw��A# '���	!3Ͼ�FB��Z���.���$e��*r,t?�����1�;}��=�4J����g��Ws�������-C�;��.�:�7�bP��ԉX���g��Z5���&�;�m���wo���{w�f�7ƍ9��v�v%G��yRxz�y.p�/f��I7B�T�q���IQ��,T����٭⒈nHY���LD^	?N)_w�O4�{�@W:�T���̈́�_iȄ��������@4p�m�y�`��
�p��gn\��M��>V�p���]��4ZƘ�r��B�ļ���ױ�avf���h��<�$��O���ٲ�hz)a���f�H賜:x��o6�����@|1Na4��W�Aa�	y:����Uw3_h���s��oHL������(�ھ=�Wl�H�_ e����x`k���p����*�(�w-v�49���*���(�"���h- ]���Z�ge�l����/$�g���"����a�W�r�E�	����%�C�����u��W�ŕj�[_�vlGy�*R�M��"�k���ީ'Ҏ 0r36�0�������c�ڟ��@����8H�u�-�}K<���K��������AU�i���9��vF��!�z�|��`4D��"c��Uc�^�'!	���b�o�2��kB׎I��x,�e�d�|�j~��f���P���{X�_���S( o���9�E��GsR� ���9���O���K����3�Q��Xx��^�DU����Z��܍&Jil1�ܓ��� n�h�Z!ʴk�KA�2\�.P>Z»�']�{�+�fMLo�Z�ڀ���bp���;Q�?m�;����Y�!-�?�(��?�5f����pk?��\�7��%Ç�%r����+E��g{��"H�7h��T��t��!�K�M�����u�%Y��I�����LO��(ٹ,�m������&�yX���;�'����m���n]��(�ұϿ�=|ì��"��x�2������M�V����%�C�����u��W�ŕj�[_�vlGy�*R�M��"�k���ީ'Ҏ 0r36�0�������c�ڟ��@����8H�u�-�}K<���K��������AU�i���9��vF��!�z�|��`4D��"c��Uc�^�'!	���b�o�2��kB׎I��x,�e�d�|�j~��f���P���{X�_���S( o���9�E��GsR� ���9���O���K����3�Q��Xx��^�DU����Z��܍&Jil1�ܓ��� n�h�Z!ʴk�KA�2\�.P>Z»�']�{�+�fMLo�Z�ڀ���bp���;Q�?m�;����Y�!-�?�(��?�5f����pk?��\�7��%Ç�%r����+E��g{��"H�7h��T��t��!�K�M�����u�%Y��I�����LO��(ٹ,�m������&�yX���;�'����m���n]��(�ұϿ�=|ì��"��x�2������M�Vʾ��	���z^���<�ˉ.�G:RL���͐^�LZw�������T)��!N�,7ؤ��U������@��+��Pb��{�##�R�ǥ�Y������]Y�:%�L2�`�Gָ����׺*�v���5�c.o�Ք�������Ve�E� r�.U�H�N�GL�%�ͤ*��tؐ���b3��������b�-�Tu婢������g}�����~�L[��]g��?�w�fݚU�#k�ϵ��C����g�'>P�Bɧ�
'I���e�8}�z`a,���俺V�~� �A����H��i:s��&�3K��9�M�=W�"�!���R�K����k|Sz��,A7��Pissl.�&aH�TH:黭L�F2c:��0^�� ��n!�0�w���1����RG�k� :��D��㚹�}[Kڪ����u��j�0�Hy��]/�tcgo7��������BS�r�����?l�@�*����5�!ԝ��>z��ɬD�0�#�ʾ��	���z^���<�ˉ.�G:RL���͐^�LZw�������T)��!N�,7ؤ��U������@��+��Pb��{�##�R�ǥ�Y������]Y�:%�L2�`�Gָ����׺*�v���5�c.o�Ք�������Ve�E� r�.U�H�N�GL�%�ͤ*��tؐ���b3��������b�-�Tu婢������g}�����~�L[��]g��?�w�fݚU�#k�ϵ��C����g�'>P�Bɧ�
'I���e�8}�z`a,���俺V�~� �A����H��i:s��&�3K��9�M�=W�"�!���R�K����k|Sz��,A7��Pissl.�&aH�TH:黭L�F2c:��0^�� ��n!�0�w���1����RG�k� :��D��㚹�}[Kڪ����u��j�0�Hy��]/�tcgo7��������BS�r�����?l�@�*����5�!ԝ��>z��ɬD�0�#���M�K)6�rNUW�*�Y��q����OI��w�����1��u6��S{�i֭Ye<�f�\[��8l�+{vy�5/�|��p�D�YQ_��1�X���V�2�I��k����y�g��ٌ$���3�xqh�9�d*���ڴ�s%.22H/�������7��*1�t� ��67��Bm�$NŁ�����l$3�
!ᕙr��B3�x�Q�9J�X/���,uT���S9Kt�=����ѷ���://�>p^7�#����m���N7v�v�~̌�)шZ+-4e)|Y�,`��
d�����[#B�i�{��J�8B.��o_*/��-ػ�e��R;��2?�=if�������.�lrq��$P� 3U�#�����Ƨs����Da氚@cY��	H�\��%Lwn��(c%�*�W;�/�t�^<7�e�&@��C�3�2bQڵ�@�*M���̹�&#��z���^Ta+�G%����I{�+,_L�.�.�����r���y����M�K)6�rNUW�*�Y��q����OI��w�����1��u6��S{�i֭Ye<�f�\[��8l�+{vy�5/�|��p�D�YQ_��1�X���V�2�I��k����y�g��ٌ$���3�xqh�9�d*���ڴ�s%.22H/�������7��*1�t� ��67��Bm�$NŁ�����l$3�
!ᕙr��B3�x�Q�9J�X/���,uT���S9Kt�=����ѷ���://�>p^7�#����m���N7v�v�~̌�)шZ+-4e)|Y�,`��
d�����[#B�i�{��J�8B.��o_*/��-ػ�e��R;��2?�=if�������.�lrq��$P� 3U�#�����Ƨs����Da氚@cY��	H�\��%Lwn��(c%�*�W;�/�t�^<7�e�&@��C�3�2bQڵ�@�*M���̹�&#��z���^Ta+�G%����I{�+,_L�.�.�����r���y��]֓ݑ�=4�1�y��8��&������n��M�=��q�����H�}����C�0���� ��	�a�%:��~1'̧���΂�D6��4�RXבg?ӭ0ȯ�|�	
A�Ұ�ė�d�]%�%�`�ncx-���.Z��f�T�����N��6= f0/���j�L-�t -|�K�d�66��'�
���ab��C��n�Ŋ(h�=�>85��<��9�ɹደ����������ګ���i�� O]f�Z4oΗn>���^����������� 8�t4N�y��4��%�"ͥwQ~ �&)#�P�c���w�@Ƈވ�ܙ�K���9��P*�'?;�'H�����7�m
�Yv�_��Y����u�N����?Ø�߀Ъ!��7!�m����{ӧ�Hݻ��4��1n�A�����+9�d̓��#�s��|1[�Nj�����{����H���f���<�n�DtC�E��`"�J��pJ��/x�Y������"��]֓ݑ�=4�1�y��8��&������n��M�=��q�����H�}����C�0���� ��	�a�%:��~1'̧���΂�D6��4�RXבg?ӭ0ȯ�|�	
A�Ұ�ė�d�]%�%�`�ncx-���.Z��f�T�����N��6= f0/���j�L-�t -|�K�d�66��'�
���ab��C��n�Ŋ(h�=�>85��<��9�ɹደ����������ګ���i�� O]f�Z4oΗn>���^����������� 8�t4N�y��4��%�"ͥwQ~ �&)#�P�c���w�@Ƈވ�ܙ�K���9��P*�'?;�'H�����7�m
�Yv�_��Y����u�N����?Ø�߀Ъ!��7!�m����{ӧ�Hݻ��4��1n�A�����+9�d̓��#�s��|1[�Nj�����{����H���f���<�n�DtC�E��`"�J��pJ��/x�Y������"����Rݤ�2�䔛؄�'��T���5��3�?F��!'�|���ΖE�K	�� 6��@B����3`�,�9�MM�ۈ�q
�Qｾ�
{�HX���T����B���~C@b2���P̭F'�@!��홿b�Fz��(��Dϴ~�k [�g���ȏT	�@i�k�K��a\,U���F��!,�Dk0���j�8+�f�u<~!9(>���m��+��_��-*L���M��א`v�^"�+���u�SpK��n�`�mo��[<^�/.�E�/��/�6�{�
-