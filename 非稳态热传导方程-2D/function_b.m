function result=function_b(x,y,t)
%%%By ������ 2021/6/9
%%%���������ڴ���Dirichlet�߽�����
if x==0
    result=exp(y+t);
elseif x==2
    result=exp(2+y+t);
elseif y==0
    result=exp(x+t);
elseif y==1
    result=exp(x+1+t);
else
    warning='û�������ı߽�����'
end