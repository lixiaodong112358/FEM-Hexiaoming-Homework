function result=function_b(x,y,t)
%%%By 李晓东 2021/6/9
%%%本函数用于处理Dirichlet边界条件
if x==0
    result=exp(y+t);
elseif x==2
    result=exp(2+y+t);
elseif y==0
    result=exp(x+t);
elseif y==1
    result=exp(x+1+t);
else
    warning='没有这样的边界条件'
end