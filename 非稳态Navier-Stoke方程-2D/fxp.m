function result=fxp(x,y,t,s)
%%%%压力的解析解
if s(1)==0&&s(2)==0
result=-(2-pi*sin(pi*x)).*cos(2*pi*y).*cos(2*pi*t);
elseif s(1)==1&&s(2)==0
    result=pi^2*cos(pi*x).*cos(2*pi*y).*cos(2*pi*t);
elseif s(1)==0&&s(2)==1
    result=-2*pi*(-2+pi*sin(pi*x)).*sin(2*pi*y).*cos(2*pi*t);
else
    warning='解析解p的导数阶输入错误，请取值0或者1'
end

end