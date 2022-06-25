function [FinalListing] = writeNeighborhoodHTMLfile(File1,NTList1,File2,NTList2,AlignedNTList1,AlignedNTList2,FileName,ErrorMsg)
HTMLtabletext='';
if length(AlignedNTList1) > 4
   [Discreps, GoodDiscreps, neighborhood1, neighborhood2] = rFindAlignmentDiscrepancies(File1,AlignedNTList1,File2,AlignedNTList2,'nearest4');
   %for html neighborhood view
   for i=1:length(neighborhood1)
        u1 = neighborhood1{i,1};
        u2 = neighborhood1{i,2};
        u3 = neighborhood1{i,3};
        u4 = neighborhood1{i,4};
        u5 = neighborhood1{i,5};
        row1 = sprintf('<tr><td>%d</td><td><label><input type="checkbox" id="index_%d_file_%d" class="jmolInline" data-coord="%s,%s,%s,%s,%s">&nbsp;</label></td><td>%s</td><td>discrepancy %8.4f</td></tr>\n',i,i,1,u1,u2,u3,u4,u5,u1, Discreps(i));
        u1 = neighborhood2{i,1};
        u2 = neighborhood2{i,2};
        u3 = neighborhood2{i,3};
        u4 = neighborhood2{i,4};
        u5 = neighborhood2{i,5};
        row2 = sprintf('<tr><td>%d</td><td><label><input type="checkbox" id="index_%d_file_%d" class="jmolInline" data-coord="%s,%s,%s,%s,%s">&nbsp;</label></td><td>%s</td></tr>\n',i,i,2,u1,u2,u3,u4,u5,u1);   
        HTMLtabletext = append(HTMLtabletext,row1,row2);
   end
   HTMLfiletemplate = fileread('r3d_align_template.html');
   HTMLfiletext = strrep(HTMLfiletemplate, '###neighborhood_list###', HTMLtabletext);
   
   fid = fopen(FileName, 'w');
   fprintf(fid, '%s', HTMLfiletext);
   fclose(fid);
elseif length(AlignedNTList1) == 4
   Discreps = xDiscrepancy(File1,AlignedNTList1,File2,AlignedNTList2);
end

