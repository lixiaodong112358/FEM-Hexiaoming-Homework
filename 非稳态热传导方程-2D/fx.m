function result=fx(x,y,t,s)
if s(1)==0&&s(2)==0
result=exp(x+y+t);
elseif s(1)==1&&s(2)==0
    result=exp(x+y+t);
elseif s(1)==0&&s(2)==1
    result=exp(x+y+t);
else
    warning='������ĵ��������������ȡֵ0����1'
end

end