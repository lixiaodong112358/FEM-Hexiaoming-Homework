function [Pb,Tb]=generate_PbTb_2D(left,bottom,h1,h2,N1,N2,P,T,basis_type)
    if basis_type==201
        Pb=P;
        Tb=T;
    elseif basis_type==202

        K=(2*N1+1)*(2*N2+1);
        Pb=zeros(2,K);
        RR=2*N2+1;  %%%每一列的行数
        Tb=zeros(6,2*N1*N2);
    %%  矩阵P
        for ii=1:K
            if mod(ii,RR)==0
                c=ii/(RR);
                r=ii-(c-1)*(RR);
            else
                r=mod(ii,RR);
                c=(ii-r)/(RR)+1;
            end

           xx=left+(c-1)*h1/2;
           yy=bottom+(r-1)*h2/2;
           Pb(1,ii)=xx;
           Pb(2,ii)=yy;
        end
    %% 矩阵 T
    k=1;
        for ce=1:N1
           for re=1:N2
               cn=ce;
               rn=re;
               j=(cn-1)*RR*2+2*rn-1;
               Tb(:,k)=[j;j+RR*2;j+2;j+RR;j+RR+1;j+1];
               k=k+1;
               Tb(:,k)=[j+2;j+RR*2;j+RR*2+2;j+RR+1;j+RR*2+1;j+RR+2];
               k=k+1;
           end
        end
    else
        warning='basis_type输入错误'
    end
end