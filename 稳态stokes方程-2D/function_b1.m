function result=function_b1(x,y)
%%%By ������ 2021/7/12
%%%���������ڴ���Dirichlet�߽�����
if x==0
    
    result=exp(-y);
    
elseif x==1
    
    result=y^2+exp(-y);
    
elseif y==-0.25
    
    result=1/16*x^2+exp(0.25);
    
elseif y==0
    
    result=1;
else
    warning='û��������Dirichlet�߽�������'
end


end