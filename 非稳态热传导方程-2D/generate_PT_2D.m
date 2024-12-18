function [P,T]=generate_PT_2D(left,bottom,h1,h2,N1,N2)
%�ú�����һ����������ֳ�N1*N2�ݣ�����N1��x�������������N2��y�����������
%�ú����γɵ�������������
%�ú�������ͨ����2021/5/4��������
K=(N1+1)*(N2+1);
P=zeros(2,K);
T=zeros(3,2*N1*N2);
%%  ����P
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
%% ���� T
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