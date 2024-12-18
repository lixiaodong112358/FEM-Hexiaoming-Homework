function result=generate_boundaryedges(N1,N2,Tb,BJ)
%%%boundaryedges
%%本函数不需要区别201，202
% clc
% clear
% % close
% N1=2;
% N2=2;
% 产生矩形区域的边界矩阵,2021/5/4,李晓东
BD=zeros(4,2*(N1+N2));
if BJ=='D'
    BD(1,:)=BD(1,:)-1;
end
BD(1,1:N1-1)=BD(1,1:N1-1)-1;    %%%底部，Neumann边界条件
% BD(1,1:end-1)=BD(1,1:end-1)-2    ;    %%%底部，Robin边界条件
%%
for i=1:N1
    BD(2,i)=1+2*N2*(i-1);
    BD(3:4,i)=Tb(1:2,BD(2,i));
end
%%
for i=(N1+1):(N1+N2)
    BD(2,i)=BD(2,N1)+1+2*(i-(N1+1));
    BD(3:4,i)=Tb(2:3,BD(2,i));
end
%%
for i=(N1+N2+1):(2*N1+N2)
    BD(2,i)=BD(2,N1+N2)-2*N2*(i-(N1+N2+1));
    BD(3:4,i)=Tb([3,1],BD(2,i));
end
%%
for i=(2*N1+N2+1):2*(N1+N2)
    BD(2,i)=BD(2,2*N1+N2)-1-2*(i-(2*N1+N2+1));
    BD(3:4,i)=Tb([3,1],BD(2,i));
end

result=BD;
end
