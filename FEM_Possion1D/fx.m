function y=fx(x,s)
if s==0
y=x*cos(x);
elseif s==1
    y=cos(x)-x*sin(x);
else
    warning='������ĵ��������������ȡֵ0����1'
end

end