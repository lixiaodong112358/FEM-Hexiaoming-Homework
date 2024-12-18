function result=function_b2(x,y,t)
%%%By 李晓东 2021/6/2
%%%本函数用于处理Dirichlet边界条件
if x==0
    result=2*cos(2*pi*t);
elseif x==1
    result=(-2/3*y^3+2)*cos(2*pi*t);
    
elseif y==-0.25
    
    result=(1/96*x+2-pi*sin(pi*x))*cos(2*pi*t);
    
elseif y==0
    
    result=(2-pi*sin(pi*x))*cos(2*pi*t);
else
    warning='没有这样的Dirichlet边界条件！'
    
end