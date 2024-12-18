function result=function_b1(x,y,t)
%%%By 李晓东 2021/7/12
%%%本函数用于处理Dirichlet边界条件
if x==0
    
    result=exp(-y)*cos(2*pi*t);
    
elseif x==1
    
    result=(y^2+exp(-y))*cos(2*pi*t);
    
elseif y==-0.25
    
    result=(1/16*x^2+exp(0.25))*cos(2*pi*t);
    
elseif y==0
    
    result=cos(2*pi*t);
else
    warning='没有这样的Dirichlet边界条件！'
end


end