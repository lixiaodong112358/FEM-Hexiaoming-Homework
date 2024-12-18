function result=fx2(x,y,t,s)
if s(1)==0&&s(2)==0
result=x.*(x-1).*y.*(y-1).*cos(t);
elseif s(1)==1&&s(2)==0
    result=(2*x-1).*(y-1).*y.*cos(t);
elseif s(1)==0&&s(2)==1
    result=(-1+x).*x.*(-1+2*y).*cos(t);
else
    warning='解析解的导数阶输入错误，请取值0或者1'
end

end