% rBarDiagram plots the alignment of File1 and File2
%   NTList1 and NTList2 are the indices of the nucleotides that are aligned

function [] = rBarDiagram(File1,NTList1,File2,NTList2,AlignedNTList1,AlignedNTList2)

View(8) = 1;                              % don't put gaps in bar diagram
View(9) = 1;                              % display the numbers

if ischar(File1),
  Filename = File1;
  File1 = zGetNTData(Filename,0);
end
if ischar(File2),
  Filename = File2;
  File2 = zGetNTData(Filename,0);
end

numBars=1;
Thickness = 0.2;                          
clf

r1=0;
r2=-.05;
for p=1:numBars
   r1 = r1-.075;                                  
   r2 = r2-.075;
   if p==1   
      A = rNumberBarDiagram(File1,NTList1,View,Thickness,r1,'above');
      text(-.5, r1+.03, File1.Filename, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 11);
   else
      A = rNumberBarDiagram(File1,NTList1,View,Thickness,r1,'none');
   end 
   
   if p==numBars
      B = rNumberBarDiagram(File2,NTList2,View,Thickness,r2,'below');
      text(-.5, r2-.035, File2.Filename, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 11);
   else
      B = rNumberBarDiagram(File2,NTList2,View,Thickness,r2,'none');
   end
     
   c = rFindAlignmentDiscrepancies(File1,AlignedNTList1,File2,AlignedNTList2,'nearest4');
   [s,t] = size(c);
   colo = c;

   if s == 1,
      for j = 1:length(i1),
         colo(j,1) = c;
      end
   end

   if (t == 1) && strcmp(class(c),'double'),
      colormap('default')
      map = colormap;
      map = map(1:55,:);
      colormap(map)
      c(c>1.1)=1.1;
      colo = map(floor(50*c),:);
   end
   i1=AlignedNTList1;
   i2=AlignedNTList2;
   for j = 1:length(i1),
      plot([A(NTList1==i1(j)) B(NTList2==i2(j))], [r1 r2],'Color',colo(j,:),'LineWidth',Thickness);
      hold on
   end
   plot([0 20],[r1 r1], 'k');
   hold on
   plot([0 20],[r2 r2], 'k');
   hold on
   text(-.5, (r1+r2)/2, 'R3D Align', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 8);
end
   
axis([-.15 20 -.1-numBars*0.075 .7-numBars*0.075]);
colorbar('southoutside','XTickLabel',{'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','>=1.1'})
axis off
orient landscape
date = regexprep(datestr(now),':', '-');
BarName=['R3D Align ' File1.Filename ' ' File2.Filename ' ' date(1:17) '.pdf'];
saveas(gcf,BarName,'pdf');