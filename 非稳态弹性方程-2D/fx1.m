function result=fx1(x,y,t,s)
if s(1)==0&&s(2)==0
result=sin(pi*x).*sin(pi*y).*sin(pi*t);
elseif s(1)==1&&s(2)==0
    result=pi*sin(pi*t).*cos(pi*x).*sin(pi*y);
elseif s(1)==0&&s(2)==1
    result=pi*sin(pi*t).*cos(pi*y).*sin(pi*x);
else
    warning='解析解的导数阶输入错误，请取值0或者1'
end

end