function x=zhuigan(A,d)
%%%%实现基于追赶法求解线性方程组
%%%%系数矩阵A必须是三对角矩阵，否则不能求解
%由于该求解方法不需要形成系数矩阵，因此计算量较低；
%数组a表示对角线下面的那一组数据，数组b表示对角线上的数据，数组c代表对角线上面的那一组数据；
%数组a的首项为0，数组b的尾项为0
a=[0;diag(A,-1)];
b=[diag(A,1);0];
c=diag(A,0);
n=max(size(c));
U=zeros(n,1);
L=zeros(n,1);
%% LU 该解法不需要形成矩阵
U(1)=b(1);
for i=2:n
    L(i)=a(i)/U(i-1);
    U(i)=b(i)-L(i)*c(i-1);
end
%% Solve y
y=zeros(n,1);
y(1)=d(1);
for i=2:n
    y(i)=d(i)-L(i)*y(i-1);
end
%% Solve x
x0=zeros(n,1);
x0(n)=y(n)/U(n);
for i=flip(1:n-1)
    x0(i)=(y(i)-c(i)*x0(i+1))/U(i);
end
x0(1)=d(1);
x0(end)=d(end);
x=x0;%输出最后的结果
end