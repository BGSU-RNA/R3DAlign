function [SM] = rMakeEdgeMatrix(VMI,List)

BI1=zeros(500000,1);
BI2=zeros(500000,1);
[m n]=size(List);
Len=zeros(m,n,'uint16');
for i = 1:m
   for j = 1:n
      Len(i,j)=length(List{i,j});
   end
end

[i j] = find(Len);

ct=0;
ct1=0;

for k=1:length(i)-1  
   b = find(logical(i(k+1:end)<=i(k)) & logical(j(k+1:end)>j(k)));
   c = find(logical(i(k+1:end)>i(k)) & logical(j(k+1:end)<=j(k)));
   d = union(b,c);
   
   for p=1:length(d)
      L1=List{i(k),j(k)};
      L2=List{i(d(p)+k),j(d(p)+k)};
      ct1=ct1+2*length(L1)*length(L2);
      for r = 1:length(L1)
         for s = 1:length(L2)
            ct=ct+1;
            len=length(BI1);
            if ct>len
               BI1(2*ct)=0;
               BI2(2*ct)=0;
            end
            BI1(ct)=L1(r);
            BI2(ct)=L2(s);
         end
      end               
   end
end

clear Len
clear i
clear j
clear List

BI1=BI1(1:ct);
BI2=BI2(1:ct);
m=length(VMI(:,1));
if isempty(BI1)
    SM=sparse(zeros(m,m));
else
    SM=sparse(BI1,BI2,ones(length(BI1),1),length(VMI),length(VMI));
end
SM=SM+SM';
SM(SM>1)=1;