function result=function_b2(x,y,t)
%%%By ������ 2021/6/2
%%%���������ڴ���Dirichlet�߽�����
if x==0
    result=2*cos(2*pi*t);
elseif x==1
    result=(-2/3*y^3+2)*cos(2*pi*t);
    
elseif y==-0.25
    
    result=(1/96*x+2-pi*sin(pi*x))*cos(2*pi*t);
    
elseif y==0
    
    result=(2-pi*sin(pi*x))*cos(2*pi*t);
else
    warning='û��������Dirichlet�߽�������'
    
end