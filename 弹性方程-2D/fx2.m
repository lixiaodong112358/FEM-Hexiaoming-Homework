function result=fx2(x,y,s)
if s(1)==0&&s(2)==0
result=x.*(x-1).*y.*(y-1);
elseif s(1)==1&&s(2)==0
    result=(-1+x).*(-1+y).*y+x.*(-1+y).*y;
elseif s(1)==0&&s(2)==1
    result=(-1+x).*x.*(-1+y)+(-1+x).*x.*y;
else
    warning='������ĵ��������������ȡֵ0����1'
end

end