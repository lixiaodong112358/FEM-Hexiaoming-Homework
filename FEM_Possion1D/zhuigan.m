function x=zhuigan(A,d)
%%%%ʵ�ֻ���׷�Ϸ�������Է�����
%%%%ϵ������A���������ԽǾ��󣬷��������
%���ڸ���ⷽ������Ҫ�γ�ϵ��������˼������ϵͣ�
%����a��ʾ�Խ����������һ�����ݣ�����b��ʾ�Խ����ϵ����ݣ�����c����Խ����������һ�����ݣ�
%����a������Ϊ0������b��β��Ϊ0
a=[0;diag(A,-1)];
b=[diag(A,1);0];
c=diag(A,0);
n=max(size(c));
U=zeros(n,1);
L=zeros(n,1);
%% LU �ýⷨ����Ҫ�γɾ���
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
x=x0;%������Ľ��
end