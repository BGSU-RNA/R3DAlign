% This file has been adopted from EvalConsensus.m and
% zCompare16SAlignments;
% Evaluate the number of near basepairs aligned with true basepairs, near with near, 
% true with true, etc.  Evaluate # of conserved nested cWW, non-nested,
% etc.  Evaluate # of nt's aligned, # of base matches.
function [] = rSetterR3DCompareSheet(MatFilename1,MatFilename2,Filename)

if exist(fullfile(pwd, 'R3D Align Output',[Filename '.xlsx'])) == 2 %#ok<EXIST>
   SpreadsheetName = fullfile(pwd, 'R3D Align Output',[Filename '.xlsx']);
elseif exist(fullfile(pwd, 'R3D Align Output',[Filename '.xls'])) == 2 %#ok<EXIST>
   SpreadsheetName = fullfile(pwd, 'R3D Align Output',[Filename '.xls']);
else
   T=cell(1,22);
   T{1,1}='File1';
   T{1,2}='File2';
   T{1,3}='Len';
   T{1,4}='Agree';
   T{1,5}='Max';
   T{1,6}='5';
   T{1,7}='10';
   for i=1:15
      T{1,i+7}=num2str(i*20);
   end
   if ispc
       fullfile(pwd, 'R3D Align Output',Filename,'.xlsx')
      xlswrite(fullfile(pwd, 'R3D Align Output',[Filename '.xlsx']),T,'Sheet1','A1');
   end
   SpreadsheetName = fullfile(pwd, 'R3D Align Output',[Filename '.xlsx']);
end
   
clear a;  clear b;
[a b]=xlsread(SpreadsheetName); %#ok<ASGLU>
T=cell(1,22);
T{1,1}=MatFilename1(1:7);
T{1,2}=MatFilename2(13:19);
L=length(b(:,1));
for i = 1:L
   if isequal(b{i,1},T{1,1}) && isequal(b{i,2},T{1,2})
      return;
   end
end
R=zCompareAlignmentVectors({MatFilename1,MatFilename2});
T{1,3}=length(R);
T{1,4}=length(find(R==0));
T{1,5}=max(abs(R));
T{1,6}=length(find(abs(R)>=5));
T{1,7}=length(find(abs(R)>=10));
for i=1:15
   T{1,i+7}=length(find(abs(R)>=i*20));
end
if ispc
   xlswrite(fullfile(pwd, 'R3D Align Output',[Filename '.xlsx']),T,'Sheet1',['A' num2str(L+1)]);
end