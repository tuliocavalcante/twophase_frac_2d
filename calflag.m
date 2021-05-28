function nflag = calflag
% determinar o flag do nó interior e fronteira de Neumann
global  coord bedge bcflag wells elem

nflag=5000*ones(size(coord,1),2);
bcflag(bcflag(:,1)<200,2)=bcflag(bcflag(:,1)<200,2)*1e-4;
bcflag(bcflag(:,1)>200,2)=bcflag(bcflag(:,1)>200,2)*(1/24)*(1/3600);

for b=1:size(bedge,1)
%     if bedge(b,1)==2 && bedge(b,2)==2
%         bedge(b,4)=101;
%         bedge(b,5)=101;
%     elseif bedge(b,1)==3 && bedge(b,2)==3
%         bedge(b,4)=102;
%         bedge(b,5)=102;
%     end
%     if bedge(b,1)==2
%         bedge(b,4)=101;
%     end
%     if bedge(b,1)==3
%         bedge(b,4)=102;
%     end
    f1=find(bcflag(:,1)==bedge(b,4));
    nflag(bedge(b,1),1)=bcflag(f1,1);
    nflag(bedge(b,1),2)=bcflag(f1,2);
    f2=find(bcflag(:,1)==bedge(b,5));
    nflag(bedge(b,2),1)=bcflag(f2,1);
    nflag(bedge(b,2),2)=bcflag(f2,2);
end

if wells(1,1)==0
    wells = double.empty(0,0);
end
    
% for i=1:size(wells,1)
%     if wells(i,5)>400 && wells(i,5)<1000
%        n1=elem(wells(i,1),1);
%        n2=elem(wells(i,1),2);
%        n3=elem(wells(i,1),3);
%        if coord(n1,1)==wells(i,7) && coord(n1,2)==wells(i,8)
%            nflag(n1,1)=190+wells(i,2);
%            nflag(n1,2)=wells(i,6);
%        elseif coord(n2,1)==wells(i,7) && coord(n2,2)==wells(i,8)
%            nflag(n2,1)=190+wells(i,2);
%            nflag(n2,2)=wells(i,6);
%        elseif coord(n3,1)==wells(i,7) && coord(n3,2)==wells(i,8)
%            nflag(n3,1)=190+wells(i,2);
%            nflag(n3,2)=wells(i,6);
%        end
%     end
% end

end