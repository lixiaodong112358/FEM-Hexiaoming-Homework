function result=function_b(x,y)
if x==-1
    result=-1.5*y*(1-y)*exp(-1+y);
elseif x==1
    result=0.5*y*(1-y)*exp(1+y);
elseif y==-1
    result=-2*x*(1-x/2)*exp(x-1);
elseif y==1
    result=0;
else
    warning='没有这样的边界条件'

end