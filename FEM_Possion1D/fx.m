function y=fx(x,s)
if s==0
y=x*cos(x);
elseif s==1
    y=cos(x)-x*sin(x);
else
    warning='解析解的导数阶输入错误，请取值0或者1'
end

end