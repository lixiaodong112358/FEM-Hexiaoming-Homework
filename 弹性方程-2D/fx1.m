function result=fx1(x,y,s)
if s(1)==0&&s(2)==0
result=sin(pi*x).*sin(pi*y);
elseif s(1)==1&&s(2)==0
    result=pi*cos(pi*x).*sin(pi*y);
elseif s(1)==0&&s(2)==1
    result=pi*cos(pi*y).*sin(pi*x);
else
    warning='������ĵ��������������ȡֵ0����1'
end

end