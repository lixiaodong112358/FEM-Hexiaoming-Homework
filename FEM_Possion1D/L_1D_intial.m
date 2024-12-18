clc
close
clear
%%%%%%初始的一维有限元程序：杨晓明 密苏里科技大学
%%%%%%Code by Xiaodong Li
%%%%%%2021/4/14
%%%%%%问题描述：-\frac{d}{d x}\left(e^{x} \frac{d u(x)}{d x}\right)=-e^{x}[\cos (x)-2 \sin (x)-x \cos (x)-x \sin (x)]
%%%%%%（0<=x<=1）
%%%%%%u(0)=0,u(1)=cos(1)
%%%%%%u=x cos(x)
%%
a=0;
b=1;
N=20;
h=(b-a)/N;
cc=@(x)exp(x);
xx=a:h:b;
%% FEM
A=zeros(N+1,N+1);
for j=1:N+1
    if j<=N
        A(j+1,j)=-1/h^2*(integral(cc,xx(j),xx(j+1)));
    end
    if j>=2
        A(j-1,j)=-1/h^2*(integral(cc,xx(j-1),xx(j)));
    end
    if j<=N&&2<=j
       A(j,j)=1/h^2*(integral(cc,xx(j-1),xx(j)))+1/h^2*(integral(cc,xx(j),xx(j+1)));
    end
end
A(1,1)=1/h^2*(integral(cc,xx(1),xx(2)));
A(N+1,N+1)=1/h^2*(integral(cc,xx(N),xx(N+1)));
%%
BB=zeros(N+1,1);
ff1=@(x,xi,h)-exp(x).*(cos(x)-2.*sin(x)-x.*cos(x)-x.*sin(x)).*(x-xi)/h;
ff2=@(x,xi,h)(-exp(x).*(cos(x)-2.*sin(x)-x.*cos(x)-x.*sin(x))).*(xi-x)/h;
for i=2:N
    %%%%%右边的积分用的是原函数表达式，不知道为什么，软件计算出现了误差；
      BB(i)=integral(@(x)ff1(x,xx(i-1),h),xx(i-1),xx(i))+ff(xx(i+1),xx(i+1),h)-ff(xx(i),xx(i+1),h);
end
% BB(1)=integral(@(x)ff2(x,x(2),h),x(1),x(2));
% BB(N+1)=integral(@(x)ff1(x,x(N),h),x(N),x(N+1));
%% Dirichlet boundary condition
A(1,:)=0;
A(1,1)=1;
A(N+1,:)=0;
A(N+1,N+1)=1;
BB(1)=0;
BB(N+1)=cos(1);
%%
yy=linsolve(full(A),BB);
y=xx.*cos(xx);
%%
plot(xx,yy,'.',xx,y,'o')