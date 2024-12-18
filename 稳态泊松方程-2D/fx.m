function result=fx(x,y,s)
if s(1)==0&&s(2)==0
result=x.*y.*(1-x/2).*(1-y).*exp(x+y);
elseif s(1)==1&&s(2)==0
    result=exp(x+y)*(1-x/2)*(1-y)*y-1/2*exp(x+y)*x*(1-y)*y+exp(x+y)*(1-x/2)*x*(1-y)*y;
elseif s(1)==0&&s(2)==1
    result=exp(x+y)*(1-x/2)*x*(1-y)-exp(x+y)*(1-x/2)*x*y+exp(x+y)*(1-x/2)*x*(1-y)*y;
else
    warning='解析解的导数阶输入错误，请取值0或者1'
end

end