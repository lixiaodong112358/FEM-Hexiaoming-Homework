function result=fxp(x,y,t,s)
%%%%ѹ���Ľ�����
if s(1)==0&&s(2)==0
result=-(2-pi*sin(pi*x)).*cos(2*pi*y).*cos(2*pi*t);
elseif s(1)==1&&s(2)==0
    result=pi^2*cos(pi*x).*cos(2*pi*y).*cos(2*pi*t);
elseif s(1)==0&&s(2)==1
    result=-2*pi*(-2+pi*sin(pi*x)).*sin(2*pi*y).*cos(2*pi*t);
else
    warning='������p�ĵ��������������ȡֵ0����1'
end

end