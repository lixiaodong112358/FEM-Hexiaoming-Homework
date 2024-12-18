function result=function_b2(x,y,t)
%%%By 李晓东 2021/6/2
%%%本函数用于处理Dirichlet边界条件
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