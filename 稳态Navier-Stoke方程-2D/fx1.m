function result=fx1(x,y,s)
%% 位移u1的解析解
if s(1)==0&&s(2)==0
result=x.^2.*y.^2+exp(-y);
elseif s(1)==1&&s(2)==0
    result=2*x.*y.^2;
elseif s(1)==0&&s(2)==1
    result=-exp(-y)+2*x.^2.*y;
else
    warning='解析解的导数阶输入错误，请取值0或者1'
end

end