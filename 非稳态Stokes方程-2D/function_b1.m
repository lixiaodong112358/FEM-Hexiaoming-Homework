function result=function_b1(x,y,t)
%%%By ������ 2021/7/12
%%%���������ڴ���Dirichlet�߽�����
if x==0
    
    result=exp(-y)*cos(2*pi*t);
    
elseif x==1
    
    result=(y^2+exp(-y))*cos(2*pi*t);
    
elseif y==-0.25
    
    result=(1/16*x^2+exp(0.25))*cos(2*pi*t);
    
elseif y==0
    
    result=cos(2*pi*t);
else
    warning='û��������Dirichlet�߽�������'
end


end