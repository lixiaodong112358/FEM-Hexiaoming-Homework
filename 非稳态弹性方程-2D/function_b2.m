function result=function_b2(x,y,t)
%%%By ������ 2021/6/2
%%%���������ڴ���Dirichlet�߽�����
if     x==-1
   result=2*(-1+y)*y*cos(t);
elseif x==1
    result=0;
elseif y==-1
    result=2*(-1+x)*x*cos(t);
elseif y==1
    result=0;
end
end