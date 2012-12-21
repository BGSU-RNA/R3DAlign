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
      elseif isequal(b{i,5}(̲�4@��ͱ�P��!�Sq ��2Ɍ�B�@Hxe���7���Ì;�eTt��� �K`6<KU���T�~�k�b�g�m-i���ˢ��ͮ��i�#�9t�w��Ӎ��ݰ3#`J4��JMY�_V'X����������PwZ�����&���(��K=�ۗʋ��xE��~�)���N����gO����k48������-�\�.�	;�L��H�?g-���\�d�>Q�9�&�X��p�'W�`	�ĝ[���Xɸ����K(]��Mo��	�#����X�vm,ЫJ���/s.���㯱U�
�Q	����q��
˗�%���`<�$��A�p@:ym��dwdvj�rރm'�+��$�25x�$9A�xO�v�<6~� {>_�7?���-�u��#6(i�cX��lɆ0��E�	�i�)����2���5M��u����O�t+�k>�g��D���4,)�%/YdW�}D	 X��^K�`���V�� ��}-����M�̲�4@��ͱ�P��!�Sq ��2Ɍ�B�@Hxe���7���Ì;�eTt��� �K`6<KU���T�~�k�b�g�m-i���ˢ��ͮ��i�#�9t�w��Ӎ��ݰ3#`J4��JMY�_V'X����������PwZ�����&���(��K=�ۗʋ��xE��~�)���N����gO����k48������-�\�.�	;�L��H�?g-���\�d�>Q�9�&�X��p�'W�`	�ĝ[���Xɸ����K(]��Mo��	�#����X�vm,ЫJ���/s.���㯱U�
�Q	����q��
˗�%���`<�$��A�p@:ym��dwdvj�rރm'�+��$�25x�$9A�xO�v�<6~� {>_�7?���-�u��#6(i�cX��lɆ0��E�	�i�)����2���5M��u����O�t+�k>�g��D���4,)�%/YdW�}D	 X��^K�`���V�� ��}-����M��KtƹZ|S˼&@��3����k�	�uu�j(ǐ���n�"
Z*FO �N��9�l�a�vrn�"�cyF�{����}9j��*f�u�F-!�S���͛󥛏�rD���h������<���F N6�Bf�}3�2�h��Hs�]���I���3T�X�~��u��7b8w&�#�{�i0�
e����	ү�>e�M;�m[��w��]�)tVo`Šb�������0f�7 �jj��M�wD�z#����5R��$2�ho�s&=�~�J�>�����\�&_���n�*�f��#�9��$�Y�$o��[�%ݐ�FQ))���~4�R����h�����t>���;�	c�Ґ	���=wOW�h��8�>�N����ܸ.E� ��}�������T7i��19�&6a���y;��cM%���@��z�y�I�/��q��eA��R��+���>-��g9u��!��lyS��6"b��h�{��>��;�t|%�ū�f��F�{�ģ��KtƹZ|S˼&@��3����k�	�uu�j(ǐ���n�"
Z*FO �N��9�l�a�vrn�"�cyF�{����}9j��*f�u�F-!�S���͛󥛏�rD���h������<���F N6�Bf�}3�2�h��Hs�]���I���3T�X�~��u��7b8w&�#�{�i0�
e����	ү�>e�M;�m[��w��]�)tVo`Šb�������0f�7 �jj��M�wD�z#����5R��$2�ho�s&=�~�J�>�����\�&_���n�*�f��#�9��$�Y�$o��[�%ݐ�FQ))���~4�R����h�����t>���;�	c�Ґ	���=wOW�h��8�>�N����ܸ.E� ��}�������T7i��19�&6a���y;��cM%���@��z�y�I�/��q��eA��R��+���>-��g9u��!��lyS��6"b��h�{��>��;�t|%�ū�f��F�{�ģ��#UB;P��Z��irKU>3�Q��DK#�Z � ����)�ʠ�jϟ_H��:�}Eb[=5�J���l�
*�vS��5$��o�������A���ܒd��+<j�[�-���xQ��K��������ǯB�+�$��(V-�؎�jU��!E � S�SO�A`�fl�ahO!nYE���?���HMm�q���$[`�$�x&A�&U+����9�����<��s�U��kC��,�*��h��3hE��s����OB$`�x�:Le�(�5ׄ��:��XL���t�V��
� G����(���w�çP ��'4�s0�H�7��TA��5ls���<����/g�,�/d��P}�r��
���v��M���b,�'��U�6���B�i/���e��]�|�(�w�O�4��"W�͚���� �sd�9��8vM�w��)ڔwWgQ�"CZ� Qt�Tk�(5g���~h9��n �K�oK�.=�W����#UB;P��Z��irKU>3�Q��DK#�Z � ����)�ʠ�jϟ_H��:�}Eb[=5�J���l�
*�vS��5$��o�������A���ܒd��+<j�[�-�