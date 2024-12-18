function [P,T]=generate_PT_2D(left,bottom,h1,h2,N1,N2)
%该函数将一个矩形区域分成N1*N2份，其中N1是x方向的网格数，N2是y方向的网格数
%该函数形成的是三角形网格
%该函数测试通过，2021/5/4；李晓东
K=(N1+1)*(N2+1);
P=zeros(2,K);
T=zeros(3,2*N1*N2);
%%  矩阵P
for ii=1:K
    if mod(ii,N2+1)==0
        c=ii/(N2+1);
        r=ii-(c-1)*(N2+1);
    else
        r=mod(ii,N2+1);
        c=(ii-r)/(N2+1)+1;
    end
   
   xx=left+(c-1)*h1;
   yy=bottom+(r-1)*h2;
   P(1,ii)=xx;
   P(2,ii)=yy;
end
%% 矩阵 T
k=1;
for ce=1:N1
   for re=1:N2
       cn=ce;
       rn=re;
       j=(cn-1)*(N2+1)+rn;
       T(:,k)=[j;j+N2+1;j+1];
       k=k+1;
       T(:,k)=[j+1;j+N2+1;j+N2+2];
       k=k+1;
   end
end
end